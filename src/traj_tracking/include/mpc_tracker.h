#pragma once
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include "OsqpEigen/OsqpEigen.h"

#define kEps 2.333e-10
namespace simple_ackermann {
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
  INVALID_STEER_RATE,
  INVALID_STEER_ANGLE,
};
struct TrackerParam {
  // mpc parameters
  int horizon_;
  double interval_;
  int state_size_;
  int input_size_;
  double weight_x_error_, weight_y_error_, weight_theta_error_;
  double weight_v_, weight_omega_;

  // vehicle parameters
  double max_vel_;
  double min_vel_;
  double max_acc_;
  double min_acc_;
  double steer_angle_rate_limit_;
  double min_turn_radius_;
  double track_width_;
  double wheel_base_;

  TrackerParam()
      : horizon_(20),
        interval_(0.05),
        state_size_(3),
        input_size_(2),
        weight_x_error_(2333.3),
        weight_y_error_(2333.3),
        weight_theta_error_(233.3),
        weight_v_(63.33),
        weight_omega_(63.33),
        max_vel_(1.5),
        min_vel_(-1.5),
        max_acc_(1.5),
        min_acc_(-1.5),
        steer_angle_rate_limit_(M_PI * 2),
        min_turn_radius_(0.3),
        track_width_(0.6),
        wheel_base_(1.0) {}
  TrackerParam(int horizon, double interval, int state_size, int input_size,
               double weight_x_error, double weight_y_error,
               double weight_theta_error, double weight_v, double weight_omega,
               double max_vel, double min_vel, double max_acc, double min_acc,
               double steer_angle_rate_limit, double min_turn_radius,
               double track_width, double wheel_base)
      : horizon_(horizon),
        interval_(interval),
        state_size_(state_size),
        input_size_(input_size),
        weight_x_error_(weight_x_error),
        weight_y_error_(weight_y_error),
        weight_theta_error_(weight_theta_error),
        weight_v_(weight_v),
        weight_omega_(weight_omega),
        max_vel_(max_vel),
        min_vel_(min_vel),
        max_acc_(max_acc),
        min_acc_(min_acc),
        steer_angle_rate_limit_(steer_angle_rate_limit),
        min_turn_radius_(min_turn_radius),
        track_width_(track_width),
        wheel_base_(wheel_base) {}
};
class MpcTracker {
 public:
  typedef std::unique_ptr<MpcTracker> UniquePtr;
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
  Vector3d x_lb_;
  Vector3d x_ub_;
  Vector2d u_lb_;
  Vector2d u_ub_;
  DMatrix A_steer_;
  DVector lb_steer_;
  DVector ub_steer_;
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
    u_lb_ << param_.min_vel_, -std::numeric_limits<double>::infinity();
    u_ub_ << param_.max_vel_, std::numeric_limits<double>::infinity();
  }
  void setSteerAngleConstraints() {
    A_steer_ = DMatrix::Zero(param_.horizon_ * 2, qp_state_size_);
    lb_steer_ = DVector::Zero(param_.horizon_ * 2);
    ub_steer_ = DVector::Zero(param_.horizon_ * 2);
    for (size_t i = 0; i < param_.horizon_; ++i) {
      A_steer_(2 * i, (param_.horizon_ + 1) * param_.state_size_ +
                          i * param_.input_size_) =
          1 / (param_.min_turn_radius_);
      A_steer_(2 * i, (param_.horizon_ + 1) * param_.state_size_ +
                          i * param_.input_size_ + 1) = 1;
      lb_steer_(2 * i) = 0.0;
      ub_steer_(2 * i) = std::numeric_limits<double>::infinity();

      A_steer_(2 * i + 1, (param_.horizon_ + 1) * param_.state_size_ +
                              i * param_.input_size_) =
          -1 / (param_.min_turn_radius_);
      A_steer_(2 * i + 1, (param_.horizon_ + 1) * param_.state_size_ +
                              i * param_.input_size_ + 1) = 1;
      lb_steer_(2 * i + 1) = -std::numeric_limits<double>::infinity();
      ub_steer_(2 * i + 1) = 0;
    }
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

  DMatrix steerRateConstraintsMatrixCaster();
  DMatrix steerRateConstraintsBoundCaster();
  DMatrix accelerateConstraintsMatrixCaster();
  DMatrix accelerateConstraintsBoundCaster();

  void castProblemToQpForm();
  SolveStatus posterioriCheck(const DVector &solution);

 public:
  MpcTracker(const TrackerParam &param);

  bool update(const Vector3d &init_state, const Trajectory2D &reference_traj);
  SolveStatus solve(Vector2d &solution);
  void getReferenceStateAndInputSeq(Trajectory3D &refer_state_seq,
                                    Trajectory2D &refer_input_seq);
  void getCurrentReferStateAndInput(Vector3d &refer_state,
                                    Vector2d &refer_input);
};

}  // namespace simple_ackermann