#pragma once
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include "OsqpEigen/OsqpEigen.h"

#define kEps 2.333e-33
namespace willand_ackermann {
enum class SolveStatus {
  SUCCESS = 0x00000000,
  // solver inner error
  SOLVER_INIT_ERROR = 0x00000001,
  SOLVER_INNER_ERROR,
  SOLVER_SET_HESSIAN_ERROR,
  SOLVER_SET_GRADIENT_ERROR,
  SOLVER_SET_CONSTRAINT_MATRIX_ERROR,
  SOLVER_SET_CONSTRAINT_BOUND_ERROR,
  // infeasible solution
  INVALID_SPEED = 0x00010000,
  INVALID_ACC,
  INVALID_STEER_ANGLE,
  INVALID_STEER_RATE,
};
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
  double weight_x_error_, weight_y_error_, weight_theta_error_;
  double weight_v_, weight_omega_;
  // vehicle parameters
  double track_width_;
  double dist_front_to_rear_;

  TrackerParam()
      : horizon_(20),
        interval_(0.2),
        state_size_(4),
        input_size_(2),
        speed_limit_(1.0),
        acc_limit_(1.0),
        front_wheel_angle_limit_(M_PI / 4),
        front_wheel_angle_rate_limit_(M_PI / 8),
        weight_x_error_(1.0),
        weight_y_error_(1.0),
        weight_theta_error_(1.0),
        weight_v_(1.0),
        weight_omega_(1.0),
        track_width_(0.5),
        dist_front_to_rear_(0.8) {}
  TrackerParam(int horizon, double interval, int state_size, int input_size,
               double speed_limit, double acc_limit,
               double front_wheel_angle_limit,
               double front_wheel_angle_rate_limit, double weight_x_error,
               double weight_y_error, double weight_theta_error,
               double weight_v, double weight_omega, double track_width,
               double dist_front_to_rear);
};
class TrajectoryTracker {
 public:
  typedef std::unique_ptr<TrajectoryTracker> UniquePtr;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DMatrix;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::VectorXd DVector;
  typedef Eigen::Vector2d Vector2d;
  typedef Eigen::Vector3d Vector3d;
  typedef std::vector<Eigen::Vector2d> Trajectory2D;
  typedef std::vector<Eigen::Vector3d> Trajectory3D;

 private:
  TrackerParam param_;
  Vector3d init_state_;
  DMatrix Q_, R_;  // weight matrices
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
  Vector3d x_lb_;
  Vector3d x_ub_;
  Vector2d u_lb_;
  Vector2d u_ub_;
  DMatrix A_steer_rate_;
  DVector lb_steer_rate_;
  DVector ub_steer_rate_;
  DMatrix A_accelerate_;
  DVector lb_accelerate_;
  DVector ub_accelerate_;
  // parameters of qp
  int qp_state_size_;
  SparseMatrix H_;
  DVector g_;
  SparseMatrix M_;  // constraint matrices
  DVector lb_, ub_;

  Trajectory3D refer_state_seq_;
  Trajectory2D refer_input_seq_;

 private:
  void calcOsqpHession();
  void calcOsqpGradient();
  void calcOsqpConstraintMatrix();
  void calcOsqpConstraintBound();
  bool setInitialState(const Vector3d &init_state) {
    if (init_state.rows() != param_.state_size_) {
      std::cout << "Invalid initial state" << std::endl;
      return false;
    }
    init_state_ = init_state;
    return true;
  }
  bool setReferenceTrajectory(const Trajectory2D &refer_traj);
  void setWeightMatrices() {
    Q_.resize(param_.state_size_, param_.state_size_);
    R_.resize(param_.input_size_, param_.input_size_);
    Q_ << param_.weight_x_error_, 0.0, 0.0, 0.0, param_.weight_y_error_, 0.0,
        0.0, 0.0, param_.weight_theta_error_;
    R_ << param_.weight_v_, 0.0, 0.0, param_.weight_omega_;
  }

  DMatrix dynamicStateMatrixCaster(const Vector3d &state,
                                   const Vector2d &input);
  DMatrix dynamicInputMatrixCaster(const Vector3d &state,
                                   const Vector2d &input);
  DVector dynamicVectorCaster(const Vector3d &state, const Vector2d &input);
  DMatrix steerRateConstraintsMatrixCaster();
  DMatrix steerRateConstraintsBoundCaster();
  DMatrix accelerateConstraintsMatrixCaster();
  DMatrix accelerateConstraintsBoundCaster();
  void setGeneralEqualityConstraints() {}

  void setGeneralInequalityConstraints() {
    // if angle limit larger than PI/2, then the constraints are not need
    if (param_.front_wheel_angle_limit_ - M_PI / 2 > -kEps ||
        param_.front_wheel_angle_limit_ + M_PI / 2 < kEps) {
      // std::cout << "Front wheel angle limit too large, constraints
      // abandoned"
      //           << std::endl;
      return;
    }
    A_inequal_ = DMatrix::Zero(4, param_.state_size_);
    B_inequal_ = DMatrix::Zero(4, param_.input_size_);
    B_inequal_ << -2 * std::tan(param_.front_wheel_angle_limit_),
        -(2 * param_.dist_front_to_rear_ -
          param_.track_width_ * std::tan(param_.front_wheel_angle_limit_)),
        -2 * std::tan(param_.front_wheel_angle_limit_),
        (2 * param_.dist_front_to_rear_ +
         param_.track_width_ * std::tan(param_.front_wheel_angle_limit_)),
        -2 * std::tan(param_.front_wheel_angle_limit_),
        -(2 * param_.dist_front_to_rear_ +
          param_.track_width_ * std::tan(param_.front_wheel_angle_limit_)),
        -2 * std::tan(param_.front_wheel_angle_limit_),
        (2 * param_.dist_front_to_rear_ -
         param_.track_width_ * std::tan(param_.front_wheel_angle_limit_));

    K_inequal_lb_.resize(4);
    K_inequal_lb_.setConstant(-std::numeric_limits<double>::infinity());
    K_inequal_ub_.resize(4);
    K_inequal_ub_.setZero();
  }
  void setGeneralBoundBoxConstraints() {
    x_lb_ = Vector3d::Zero();
    x_ub_ = Vector3d::Zero();
    x_lb_ << -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity();
    x_ub_ << +std::numeric_limits<double>::infinity(),
        +std::numeric_limits<double>::infinity(),
        +std::numeric_limits<double>::infinity();
    u_lb_ = Vector2d::Zero();
    u_ub_ = Vector2d::Zero();
    u_lb_ << -param_.speed_limit_, -std::numeric_limits<double>::infinity();
    u_ub_ << param_.speed_limit_, std::numeric_limits<double>::infinity();
  }
  void setSteerRateConstraints() {
    A_steer_rate_ = steerRateConstraintsMatrixCaster();
    DMatrix b = steerRateConstraintsBoundCaster();
    lb_steer_rate_ = b.col(0);
    ub_steer_rate_ = b.col(1);
  }
  void setAccelerateConstraints() {
    A_accelerate_ = accelerateConstraintsMatrixCaster();
    DMatrix b = accelerateConstraintsBoundCaster();
    lb_accelerate_ = b.col(0);
    ub_accelerate_ = b.col(1);
  }
  void castProblemToQpForm();
  SolveStatus posterioriCheck(const DVector &solution);

 public:
  TrajectoryTracker(const TrackerParam &param);

  bool update(const Vector3d &init_state, const Trajectory2D &reference_traj);
  SolveStatus solve(DVector &solution);
  void getReferenceStateAndInputSeq(Trajectory3D &refer_state_seq,
                                    Trajectory2D &refer_input_seq);
  void getCurrentReferStateAndInput(Vector3d &refer_state,
                                    Vector2d &refer_input);
};

}  // namespace willand_ackermann