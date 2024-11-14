#pragma once
#include "Eigen/Dense"
#include "OsqpEigen/OsqpEigen.h"
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <cmath>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#define kEps 2.22507e-308
namespace willand_ackermann {
struct TrackerParam {
  // mpc parameters
  int horizon_;
  double interval_;
  int state_size_;
  int input_size_;
  double speed_limit_;
  double acc_limit_;
  double front_wheel_angle_limit_;
  double front_wheel_angle_rate_limit_;
  // vehicle parameters
  double track_width_;
  double dist_front_to_rear_;

  TrackerParam()
      : horizon_(20), interval_(0.2), state_size_(4), input_size_(2),
        speed_limit_(1.0), acc_limit_(1.0), front_wheel_angle_limit_(M_PI / 4),
        front_wheel_angle_rate_limit_(M_PI / 8), track_width_(0.5),
        dist_front_to_rear_(0.8) {}
  TrackerParam(int horizon, double interval, int state_size, int input_size,
               double speed_limit, double acc_limit,
               double front_wheel_angle_limit,
               double front_wheel_angle_rate_limit, double track_width,
               double dist_front_to_rear);
};
class TrajectoryTracker {
public:
  typedef std::unique_ptr<TrajectoryTracker> UniquePtr;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DMatrix;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::VectorXd DVector;
  typedef Eigen::Vector2d Point2d;
  typedef Eigen::Vector2d Vector2d;
  typedef std::vector<Eigen::Vector2d> Trajectory2D;
  typedef std::vector<Eigen::VectorXd> TrajectoryXD;
  typedef std::function<DMatrix(const TrackerParam &param, const DVector &,
                                const DVector &)>
      MatrixCaster;
  typedef std::function<DVector(const TrackerParam &param, const DVector &,
                                const DVector &)>
      VectorCaster;
  typedef std::function<DMatrix(const TrackerParam &param, TrajectoryXD,
                                TrajectoryXD)>
      UserCustomizeMatrixCaster;
  typedef std::function<DMatrix(const TrackerParam &param, TrajectoryXD,
                                TrajectoryXD)>
      UserCustomizeBoundCaster;

private:
  TrackerParam param_;
  DVector init_state_;
  DMatrix Q_, R_; // weight matrices
  DMatrix Ad_;
  DMatrix Bd_;
  DVector Kd_;
  DMatrix A_equal_;
  DMatrix B_equal_;
  DVector K_equal_;
  DMatrix A_inequal_;
  DMatrix B_inequal_;
  DVector K_inequal_lb_;
  DVector K_inequal_ub_;
  DVector x_lb_;
  DVector x_ub_;
  DVector u_lb_;
  DVector u_ub_;
  DMatrix P_customize_;
  DVector lb_customize_;
  DVector ub_customize_;
  // parameters of qp
  int qp_state_size_;
  SparseMatrix H_;
  DVector g_;
  SparseMatrix M_; // constraint matrices
  DVector lb_, ub_;
  // DMatrix cons_bd_;
  // std::vector<DMatrix> Ad_seq_;
  // std::vector<DMatrix> Bd_seq_;
  // std::vector<DMatrix> Kd_seq_;
  // const Trajectory2D *refer_traj_ptr_;
  TrajectoryXD refer_state_seq_;
  TrajectoryXD refer_input_seq_;
  MatrixCaster DynamicStateMatrixCaster;
  MatrixCaster DynamicInputMatrixCaster;
  VectorCaster DynamicVectorCaster;

private:
  void calcOsqpHession();
  void calcOsqpGradient();
  void calcOsqpConstraintMatrix();
  void calcOsqpConstraintBound();
  inline bool setInitialState(const DVector &init_state) {
    if (init_state.rows() != param_.state_size_) {
      std::cout << "Invalid initial state" << std::endl;
      return false;
    }
    init_state_ = init_state;
    return true;
  }
  inline bool setWeightMatrices(const DMatrix &Q, const DMatrix &R) {
    if (Q.rows() != param_.state_size_ || Q.cols() != param_.state_size_ ||
        R.rows() != param_.input_size_ || R.cols() != param_.input_size_) {
      std::cout << "Invalid weight matrix" << std::endl;
      return false;
    }
    Q_ = Q;
    R_ = R;
    return true;
  }
  inline void
  setDynamicParamsCaster(const MatrixCaster &dynamic_state_matrix_caster,
                         const MatrixCaster &dynamic_input_matrix_caster,
                         const VectorCaster &dynamic_vector_caster) {
    DynamicStateMatrixCaster = dynamic_state_matrix_caster;
    DynamicInputMatrixCaster = dynamic_input_matrix_caster;
    DynamicVectorCaster = dynamic_vector_caster;
  }

  inline void setGeneralEqualityConstraints(const DMatrix &A, const DMatrix &B,
                                            const DVector &K) {
    A_equal_ = A;
    B_equal_ = B;
    K_equal_ = K;
    // return true;
  }

  inline void setGeneralInequalityConstraints(const DMatrix &A,
                                              const DMatrix &B,
                                              const DVector &lb,
                                              const DVector &ub) {
    A_inequal_ = A;
    B_inequal_ = B;
    K_inequal_lb_ = lb;
    K_inequal_ub_ = ub;
    // return true;
  }
  inline void setGeneralBoundBoxConstraints(const DVector &x_lb,
                                            const DVector &x_ub,
                                            const DVector &u_lb,
                                            const DVector &u_ub) {
    x_lb_ = x_lb;
    x_ub_ = x_ub;
    u_lb_ = u_lb;
    u_ub_ = u_ub;
  }
  bool setReferenceTrajectory(const Trajectory2D &refer_traj);
  void CastProblemToQpForm();

public:
  TrajectoryTracker(const TrackerParam &param);

  bool init(const DMatrix &Q, const DMatrix &R,
            const MatrixCaster &dynamic_state_matrix_caster,
            const MatrixCaster &dynamic_input_matrix_caster,
            const VectorCaster &dynamic_vector_caster, const DMatrix &A_equal,
            const DMatrix &B_equal, const DVector &K_equal,
            const DMatrix &A_inequal, const DMatrix &B_inequal,
            const DVector &K_inequal_lb, const DVector &K_inequal_ub,
            const DVector &x_lb, const DVector &x_ub, const DVector &u_lb,
            const DVector &u_ub, const DVector &init_state,
            const Trajectory2D &reference_traj);
  bool solve(DVector &solution);
  void
  addUserCustomizedConstraints(const UserCustomizeMatrixCaster &matrix_caster,
                               const UserCustomizeBoundCaster &bound_caster);
  void printRefereceStateSeq();
  void printRefereceInputSeq();
};

} // namespace willand_ackermann