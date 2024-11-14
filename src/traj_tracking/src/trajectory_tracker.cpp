#include "trajectory_tracker.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>

namespace willand_ackermann {

TrackerParam::TrackerParam(int horizon, double interval, int state_size,
                           int input_size, double speed_limit, double acc_limit,
                           double front_wheel_angle_limit,
                           double front_wheel_angle_rate_limit,
                           double track_width, double dist_front_to_rear)
    : horizon_(horizon), interval_(interval), state_size_(state_size),
      input_size_(input_size), speed_limit_(speed_limit), acc_limit_(acc_limit),
      front_wheel_angle_limit_(front_wheel_angle_rate_limit),
      front_wheel_angle_rate_limit_(front_wheel_angle_rate_limit),
      track_width_(track_width), dist_front_to_rear_(dist_front_to_rear) {}

TrajectoryTracker::TrajectoryTracker(const TrackerParam &param)
    : param_(param) {
  // all states and inputs are stacked into a large vector
  qp_state_size_ = param_.state_size_ * (param_.horizon_ + 1) +
                   param_.input_size_ * param_.horizon_;
  // pre-allocate matrices memory
  H_.resize(qp_state_size_, qp_state_size_);
  g_.resize(qp_state_size_);
};

bool TrajectoryTracker::init(const DMatrix &Q, const DMatrix &R,
                             const MatrixCaster &dynamic_state_matrix_caster,
                             const MatrixCaster &dynamic_input_matrix_caster,
                             const VectorCaster &dynamic_vector_caster,
                             const DMatrix &A_equal, const DMatrix &B_equal,
                             const DVector &K_equal, const DMatrix &A_inequal,
                             const DMatrix &B_inequal,
                             const DVector &K_inequal_lb,
                             const DVector &K_inequal_ub, const DVector &x_lb,
                             const DVector &x_ub, const DVector &u_lb,
                             const DVector &u_ub, const DVector &init_state,
                             const Trajectory2D &refer_traj) {

  if (!setInitialState(init_state)) {
    std::cout << "Invalid initial state!" << std::endl;
    return false;
  }
  if (!setReferenceTrajectory(refer_traj)) {
    std::cout << "Invalid reference trajectory!" << std::endl;
    return false;
  }
  if (!setWeightMatrices(Q, R)) {
    std::cout << "Invalid weight matrices trajectory!" << std::endl;
    return false;
  }
  setDynamicParamsCaster(dynamic_state_matrix_caster,
                         dynamic_input_matrix_caster, dynamic_vector_caster);
  setGeneralEqualityConstraints(A_equal, B_equal, K_equal);
  setGeneralInequalityConstraints(A_inequal, B_inequal, K_inequal_lb,
                                  K_inequal_ub);
  setGeneralBoundBoxConstraints(x_lb, x_ub, u_lb, u_ub);
  // cast everything need
  CastProblemToQpForm();
  return true;
}

bool TrajectoryTracker::solve(DVector &solution) {
  // instantiate the solver
  OsqpEigen::Solver solver;

  // settings
  // solver.settings()->setVerbosity(false);
  solver.settings()->setWarmStart(true);

  // set the initial data of the QP solver
  solver.data()->setNumberOfVariables(param_.state_size_ *
                                          (param_.horizon_ + 1) +
                                      param_.input_size_ * param_.horizon_);
  if (!solver.data()->setHessianMatrix(H_)) {
    return false;
  }
  if (!solver.data()->setGradient(g_)) {
    return false;
  }

  solver.data()->setNumberOfConstraints(M_.rows());
  if (!solver.data()->setLinearConstraintsMatrix(M_)) {

    return false;
  }
  if (!solver.data()->setLowerBound(lb_)) {
    return false;
  }
  if (!solver.data()->setUpperBound(ub_)) {
    return false;
  }

  // instantiate the solver
  if (!solver.initSolver()) {
    return false;
  }
  if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) {
    std::cout << "Osqp solver inner error" << std::endl;
    return false;
  }
  // get the controller input
  solution = solver.getSolution();

  return true;
}

void TrajectoryTracker::CastProblemToQpForm() {
  calcOsqpHession();
  calcOsqpGradient();
  calcOsqpConstraintMatrix();
  calcOsqpConstraintBound();
}

void TrajectoryTracker::calcOsqpHession() {
  // weights for state variables
  DMatrix H;
  H.resize(qp_state_size_, qp_state_size_);
  for (size_t i = 0; i <= param_.horizon_; ++i) {
    H.block(i * param_.state_size_, i * param_.state_size_, param_.state_size_,
            param_.state_size_) = Q_;
  }
  // weights for input variables
  for (size_t i = 0; i < param_.horizon_; ++i) {
    H.block(param_.state_size_ * (param_.horizon_ + 1) + i * param_.input_size_,
            param_.state_size_ * (param_.horizon_ + 1) + i * param_.input_size_,
            param_.input_size_, param_.input_size_) = R_;
  }
  H_ = H.sparseView();
}

void TrajectoryTracker::calcOsqpGradient() {
  /*
   ! Note: the standard form of quadratic programming (QP) problem is as follow
   ! 1/2 * x.transpose() * P * x + q.transpose() * x
   ! that's why the gradient needs to be multiplied by 0.5
  */
  // weights for input variables are all zero
  g_.setZero();
  // weights for state variables
  for (size_t i = 0; i <= param_.horizon_; ++i) {
    g_.segment(i * param_.state_size_, param_.state_size_) =
        -Q_ * refer_state_seq_.at(i); // has multiplied by 0.5
  }
}
void TrajectoryTracker::calcOsqpConstraintMatrix() {
  // M's rows number = initial state + dynamic model + equality cons +
  // inequality cons + bounding box
  int nums_of_initial_state = param_.state_size_,
      nums_of_dynamic = param_.state_size_ * param_.horizon_,
      nums_of_equality_cons = A_equal_.rows() * param_.horizon_,
      nums_of_inequality_cons = A_inequal_.rows() * param_.horizon_,
      nums_of_state_bounding_box = x_lb_.rows() * param_.horizon_,
      nums_of_input_bounding_box = u_lb_.rows() * param_.horizon_,
      nums_of_customize_cons = P_customize_.rows();
  // pre-allocate
  DMatrix M;
  M.resize(nums_of_initial_state + nums_of_dynamic + nums_of_equality_cons +
               nums_of_inequality_cons + nums_of_state_bounding_box +
               nums_of_input_bounding_box + nums_of_customize_cons,
           qp_state_size_);

  // initial state cons
  M.block(0, 0, param_.state_size_, param_.state_size_).setIdentity();
  // dynamic
  for (size_t i = 0; i < param_.horizon_; ++i) {
    auto &refer_state = refer_state_seq_.at(i);
    auto &refer_input = refer_input_seq_.at(i);
    M.block(nums_of_initial_state + i * param_.state_size_,
            i * param_.state_size_, param_.state_size_, param_.state_size_) =
        DynamicStateMatrixCaster(param_, refer_state, refer_input);
    M.block(nums_of_initial_state + i * param_.state_size_,
            (i + 1) * param_.state_size_, param_.state_size_,
            param_.state_size_) =
        -Eigen::MatrixXd::Identity(param_.state_size_, param_.state_size_);
    M.block(nums_of_initial_state + i * param_.state_size_,
            (param_.horizon_ + 1) * param_.state_size_ + i * param_.input_size_,
            param_.state_size_, param_.input_size_) =
        DynamicInputMatrixCaster(param_, refer_state, refer_input);
  }
  // equality cons
  if (nums_of_equality_cons != 0) {
    int block_rows = A_equal_.rows();
    int start_row_offset = nums_of_initial_state + nums_of_dynamic;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      M.block(start_row_offset + i * block_rows, i * param_.state_size_,
              block_rows, param_.state_size_) = A_equal_;
      M.block(start_row_offset + i * block_rows,
              (param_.horizon_ + 1) * param_.state_size_ +
                  i * param_.input_size_,
              block_rows, param_.input_size_) = B_equal_;
    }
  } else {
    std::cout << "No equality constraints need to be casted." << std::endl;
  }
  // inequality cons
  if (nums_of_inequality_cons != 0) {
    int block_rows = A_inequal_.rows();
    int start_row_offset =
        nums_of_initial_state + nums_of_dynamic + nums_of_equality_cons;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      M.block(start_row_offset + i * block_rows, i * param_.state_size_,
              block_rows, param_.state_size_) = A_inequal_;
      M.block(start_row_offset + i * block_rows,
              (param_.horizon_ + 1) * param_.state_size_ +
                  i * param_.input_size_,
              block_rows, param_.input_size_) = B_inequal_;
    }
  } else {
    std::cout << "No inequality constraints need to be casted." << std::endl;
  }

  // state bounding box cons
  if (nums_of_state_bounding_box != 0) {
    int block_rows = param_.state_size_;
    int start_row_offset = nums_of_initial_state + nums_of_dynamic +
                           nums_of_equality_cons + nums_of_inequality_cons;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      M.block(start_row_offset + i * block_rows, i * param_.state_size_,
              block_rows, param_.state_size_) =
          Eigen::MatrixXd::Identity(param_.state_size_, param_.state_size_);
    }
  } else {
    std::cout << "No state bounding box constraints need to be casted."
              << std::endl;
  }
  // input bounding box cons
  if (nums_of_input_bounding_box != 0) {
    int block_rows = param_.input_size_;
    int start_row_offset = nums_of_initial_state + nums_of_dynamic +
                           nums_of_equality_cons + nums_of_inequality_cons +
                           nums_of_state_bounding_box;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      M.block(start_row_offset + i * block_rows,
              (param_.horizon_ + 1) * param_.state_size_ +
                  i * param_.input_size_,
              block_rows, param_.input_size_) =
          Eigen::MatrixXd::Identity(param_.input_size_, param_.input_size_);
    }
    // std::cout << "input bounding constraints start at row: " <<
    // start_row_offset
    //           << std::endl;
  } else {
    std::cout << "No input bounding box constraints need to be casted."
              << std::endl;
  }
  // user customize customized constraints
  if (nums_of_customize_cons != 0) {
    int start_row_offset = nums_of_initial_state + nums_of_dynamic +
                           nums_of_equality_cons + nums_of_inequality_cons +
                           nums_of_state_bounding_box +
                           nums_of_input_bounding_box;
    M.block(start_row_offset, 0, P_customize_.rows(), P_customize_.cols()) =
        P_customize_;
  } else {
    std::cout << "No user customized constraints need to be casted."
              << std::endl;
  }
  M_ = M.sparseView();
}
void TrajectoryTracker::calcOsqpConstraintBound() {
  int nums_of_initial_state = param_.state_size_,
      nums_of_dynamic = param_.state_size_ * param_.horizon_,
      nums_of_equality_cons = A_equal_.rows() * param_.horizon_,
      nums_of_inequality_cons = A_inequal_.rows() * param_.horizon_,
      nums_of_state_bounding_box = x_lb_.rows() * param_.horizon_,
      nums_of_input_bounding_box = u_lb_.rows() * param_.horizon_,
      nums_of_customize_cons = P_customize_.rows();
  int nums_of_cons_rows = nums_of_initial_state + nums_of_dynamic +
                          nums_of_equality_cons + nums_of_inequality_cons +
                          nums_of_state_bounding_box +
                          nums_of_input_bounding_box + nums_of_customize_cons;
  lb_.resize(nums_of_cons_rows);
  ub_.resize(nums_of_cons_rows);

  // initial state cons
  lb_.segment(0, param_.state_size_) = init_state_;
  ub_.segment(0, param_.state_size_) = init_state_;
  // dynamic cons
  for (size_t i = 0; i < param_.horizon_; ++i) {
    auto &refer_state = refer_state_seq_.at(i);
    auto &refer_input = refer_input_seq_.at(i);
    lb_.segment(nums_of_initial_state + i * param_.state_size_,
                param_.state_size_) =
        -DynamicVectorCaster(param_, refer_state, refer_input);
    ub_.segment(nums_of_initial_state + i * param_.state_size_,
                param_.state_size_) =
        -DynamicVectorCaster(param_, refer_state, refer_input);
  }
  // equality cons
  if (nums_of_equality_cons != 0) {
    int block_rows = A_equal_.rows();
    int start_row_offset = nums_of_initial_state + nums_of_dynamic;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      lb_.segment(start_row_offset + i * block_rows, block_rows) = K_equal_;
      ub_.segment(start_row_offset + i * block_rows, block_rows) = K_equal_;
    }
  }
  // inequality cons
  if (nums_of_inequality_cons != 0) {
    int block_rows = A_inequal_.rows();
    int start_row_offset =
        nums_of_initial_state + nums_of_dynamic + nums_of_equality_cons;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      lb_.segment(start_row_offset + i * block_rows, block_rows) =
          K_inequal_lb_;
      ub_.segment(start_row_offset + i * block_rows, block_rows) =
          K_inequal_ub_;
    }
  }
  // state bounding box cons
  if (nums_of_state_bounding_box != 0) {
    int block_rows = param_.state_size_;
    int start_row_offset = nums_of_initial_state + nums_of_dynamic +
                           nums_of_equality_cons + nums_of_inequality_cons;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      lb_.segment(start_row_offset + i * block_rows, block_rows) = x_lb_;
      ub_.segment(start_row_offset + i * block_rows, block_rows) = x_ub_;
    }
  }
  // input bounding box cons
  if (nums_of_input_bounding_box != 0) {
    int block_rows = param_.input_size_;
    int start_row_offset = nums_of_initial_state + nums_of_dynamic +
                           nums_of_equality_cons + nums_of_inequality_cons +
                           nums_of_state_bounding_box;
    for (size_t i = 0; i < param_.horizon_; ++i) {
      lb_.segment(start_row_offset + i * block_rows, block_rows) = u_lb_;
      ub_.segment(start_row_offset + i * block_rows, block_rows) = u_ub_;
    }
  }
  // user customized constraints
  if (P_customize_.rows() != 0) {
    int start_row_offset = nums_of_initial_state + nums_of_dynamic +
                           nums_of_equality_cons + nums_of_inequality_cons +
                           nums_of_state_bounding_box +
                           nums_of_input_bounding_box;
    lb_.segment(start_row_offset, lb_customize_.rows()) = lb_customize_;
    ub_.segment(start_row_offset, ub_customize_.rows()) = ub_customize_;
  }
}
bool TrajectoryTracker::setReferenceTrajectory(const Trajectory2D &refer_traj) {
  if (refer_traj.size() <= 1) {
    std::cout << "Invalid reference trajectory" << std::endl;
    return false;
  }

  // given a 2-d trajectory, calculate the reference states and inputs
  // calculate reference states: (x, y, theta and velocity)
  refer_state_seq_.clear();
  refer_state_seq_.reserve(param_.horizon_ + 1);
  refer_state_seq_.resize(param_.horizon_ + 1);
  for (size_t i = 0; i < param_.horizon_; ++i) {
    auto &refer_state = refer_state_seq_.at(i);
    refer_state.resize(param_.state_size_);
    refer_state.segment(0, 2) = refer_traj.at(i);
    double delta_x = refer_traj.at(i + 1)(0) - refer_traj.at(i)(0);
    double delta_y = refer_traj.at(i + 1)(1) - refer_traj.at(i)(1);
    refer_state(2) = std::atan(delta_y / (delta_x + kEps));
    refer_state(2) = delta_x > 0 ? refer_state(2)
                                 : delta_y > 0 ? refer_state(2) + M_PI
                                               : refer_state(2) - M_PI;
    refer_state(3) =
        std::sqrt(delta_x * delta_x + delta_y * delta_y) / param_.interval_;
  }
  auto &refer_state = refer_state_seq_.back();
  refer_state.resize(param_.state_size_);
  refer_state.segment(0, 2) = refer_traj.back();
  refer_state.segment(2, 2) =
      refer_state_seq_.at(param_.horizon_ - 1).segment(2, 2);

  // calculate approximate inputs: omega and accelaration
  refer_input_seq_.clear();
  refer_input_seq_.reserve(param_.horizon_);
  refer_input_seq_.resize(param_.horizon_);
  for (size_t i = 1; i < param_.horizon_; ++i) {
    auto &refer_input = refer_input_seq_.at(i);
    double acc = (refer_state_seq_.at(i + 1)(3) - refer_state_seq_.at(i)(3)) /
                 param_.interval_;
    refer_input.resize(param_.input_size_);
    refer_input(1) = acc;
    // calculate curvature by three points, assume moving along a circle path
    // in a short distance
    Point2d A = refer_state_seq_.at(i - 1).segment<2>(0);
    Point2d B = refer_state_seq_.at(i).segment<2>(0);
    Point2d C = refer_state_seq_.at(i + 1).segment<2>(0);
    Vector2d ab = B - A, ac = C - A, bc = C - B;
    double angle_included =
        std::acos(ab.dot(ac) / (ab.norm() * ac.norm() + kEps));
    double radius = std::abs(bc.norm() / 2 / std::sin(angle_included));
    double omega = refer_state_seq_.at(i)(3) / radius;
    refer_input(0) = omega;
  }
  // the first input vector is not yet be calculated
  refer_input_seq_.front() = refer_input_seq_.at(1);
  return true;
}
void TrajectoryTracker::addUserCustomizedConstraints(
    const UserCustomizeMatrixCaster &matrix_caster,
    const UserCustomizeBoundCaster &bound_caster) {
  DMatrix m = matrix_caster(param_, refer_state_seq_, refer_input_seq_);
  DMatrix b = bound_caster(param_, refer_state_seq_, refer_input_seq_);
  if (P_customize_.rows() != 0) {
    P_customize_.conservativeResize(P_customize_.rows() + m.rows(),
                                    P_customize_.cols());
    lb_customize_.conservativeResize(lb_customize_.rows() + b.rows());
    ub_customize_.conservativeResize(ub_customize_.rows() + b.rows());
  } else {
    P_customize_.resize(m.rows(), m.cols());
    lb_customize_.resize(b.rows());
    ub_customize_.resize(b.rows());
  }
  P_customize_.bottomRows(m.rows()) = m;
  lb_customize_.tail(b.rows()) = b.col(0);
  ub_customize_.tail(b.rows()) = b.col(1);
}
void TrajectoryTracker::printRefereceStateSeq() {
  Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "", "(", ")");
  for (size_t i = 0; i <= param_.horizon_; ++i) {
    std::cout << std::fixed << std::setprecision(2) << std::setw(4)
              << std::showpos << "time stamp = " << i * param_.interval_
              << ", reference state = "
              << refer_state_seq_.at(i).transpose().format(CleanFmt)
              << std::endl;
  }
}
void TrajectoryTracker::printRefereceInputSeq() {
  Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "", "(", ")");
  for (size_t i = 0; i < param_.horizon_; ++i) {
    std::cout << std::fixed << std::setprecision(2) << std::setw(4)
              << std::showpos << "time stamp = " << i * param_.interval_
              << ", reference state = "
              << refer_input_seq_.at(i).transpose().format(CleanFmt)
              << std::endl;
  }
}

}; // namespace willand_ackermann