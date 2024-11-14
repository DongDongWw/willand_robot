#include "trajectory_tracker.h"
#include <Eigen/Dense>
#include <Eigen/src/Core/IO.h>
#include <Eigen/src/Core/Matrix.h>
#include <bits/c++config.h>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <linux/limits.h>
#include <memory>
using namespace willand_ackermann;

int main() {
  const int horizon = 40;                                   // duration = 8 secs
  const double interval = 0.2;                              // unit, sec
  const int state_size = 4;                                 // (x, y, theta, v)
  const int input_size = 2;                                 // (omega, acc)
  constexpr double speed_limit = 2.0;                       // unit, m
  constexpr double acc_limit = 2.0;                         // unit, m
  constexpr double front_wheel_angle_limit = M_PI / 4;      // unit, rad
  constexpr double front_wheel_angle_rate_limit = M_PI / 8; // unit, rad per sec
  constexpr double track_width = 0.5;                       // unit, m
  constexpr double dist_front_to_rear = 0.8;                // unit, m

  TrackerParam param(horizon, interval, state_size, input_size, speed_limit,
                     acc_limit, front_wheel_angle_limit,
                     front_wheel_angle_rate_limit, track_width,
                     dist_front_to_rear);

  TrajectoryTracker::DMatrix Q;
  TrajectoryTracker::DMatrix R;
  TrajectoryTracker::DMatrix A_equal;
  TrajectoryTracker::DMatrix B_equal;
  TrajectoryTracker::DVector K_equal;
  TrajectoryTracker::DMatrix A_inequal;
  TrajectoryTracker::DMatrix B_inequal;
  TrajectoryTracker::DVector K_inequal_lb;
  TrajectoryTracker::DVector K_inequal_ub;
  TrajectoryTracker::DVector x_lb;
  TrajectoryTracker::DVector x_ub;
  TrajectoryTracker::DVector u_lb;
  TrajectoryTracker::DVector u_ub;
  TrajectoryTracker::DVector init_state;

  Q.resize(state_size, state_size);
  R.resize(input_size, input_size);
  Q << 10, 0, 0, 0, 0, 10, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
  R << 1, 0, 0, 1;
  std::cout << Q << std::endl;
  std::cout << R << std::endl;
  A_inequal.resize(4, state_size);
  B_inequal.resize(4, input_size);
  A_inequal << 0, 0, 0, -2 * std::tan(front_wheel_angle_limit), 0, 0, 0,
      -2 * std::tan(front_wheel_angle_limit), 0, 0, 0,
      -2 * std::tan(front_wheel_angle_limit), 0, 0, 0,
      -2 * std::tan(front_wheel_angle_limit);
  B_inequal << -(2 * dist_front_to_rear -
                 track_width * std::tan(front_wheel_angle_limit)),
      0,
      (2 * dist_front_to_rear +
       track_width * std::tan(front_wheel_angle_limit)),
      0,
      -(2 * dist_front_to_rear +
        track_width * std::tan(front_wheel_angle_limit)),
      0,
      (2 * dist_front_to_rear -
       track_width * std::tan(front_wheel_angle_limit)),
      0;
  K_inequal_lb.resize(4);
  K_inequal_lb.setConstant(-std::numeric_limits<double>::infinity());
  K_inequal_ub.resize(4);
  K_inequal_ub.setZero();

  x_lb.resize(state_size);
  x_ub.resize(state_size);
  x_lb << -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(), -M_PI, -speed_limit;
  x_ub << +std::numeric_limits<double>::infinity(),
      +std::numeric_limits<double>::infinity(), +M_PI, speed_limit;
  u_lb.resize(input_size);
  u_ub.resize(input_size);
  u_lb << -std::numeric_limits<double>::infinity(), -acc_limit;
  u_ub << std::numeric_limits<double>::infinity(), acc_limit;

  init_state.resize(state_size);
  init_state << 0, 0, 0, 0;
  // funciton mapping parameters to discrete linear matrix
  auto dynamic_state_matrix_caster =
      [](const TrackerParam &param, const TrajectoryTracker::DVector &x_refer,
         const TrajectoryTracker::DVector &u_refer)
      -> TrajectoryTracker::DMatrix {
    int state_size = param.state_size_;
    int input_size = param.input_size_;
    double interval = param.interval_;
    double theta = x_refer(2), v = x_refer(3);
    TrajectoryTracker::DMatrix partial_x;
    partial_x.resize(state_size, state_size);
    partial_x << 0, 0, 0, 0, 0, 0, 0, 0, -v * std::sin(theta),
        v * std::cos(theta), 0, 0, std::cos(theta), std::sin(theta), 0, 0;
    return Eigen::MatrixXd::Identity(state_size, state_size) +
           interval * partial_x.transpose();
  };

  auto dynamic_input_matrix_caster =
      [](const TrackerParam &param, const TrajectoryTracker::DVector &x_refer,
         const TrajectoryTracker::DVector &u_refer)
      -> TrajectoryTracker::DMatrix {
    int state_size = param.state_size_;
    int input_size = param.input_size_;
    double interval = param.interval_;
    TrajectoryTracker::DMatrix partial_u;
    partial_u.resize(input_size, state_size);
    partial_u << 0, 0, 1, 0, 0, 0, 0, 1;
    return interval * partial_u.transpose();
  };

  auto dynamic_vector_caster = [](const TrackerParam &param,
                                  const TrajectoryTracker::DVector &x_refer,
                                  const TrajectoryTracker::DVector &u_refer)
      -> TrajectoryTracker::DVector {
    int state_size = param.state_size_;
    int input_size = param.input_size_;
    double interval = param.interval_;
    double theta = x_refer(2), v = x_refer(3), omege = u_refer(0),
           acc = u_refer(1);
    TrajectoryTracker::DVector x_dot(state_size);
    TrajectoryTracker::DMatrix partial_x(state_size, state_size);
    TrajectoryTracker::DMatrix partial_u(input_size, state_size);

    x_dot << v * std::cos(theta), v * std::sin(theta), omege, acc;
    partial_x << 0, 0, 0, 0, 0, 0, 0, 0, -v * std::sin(theta),
        v * std::cos(theta), 0, 0, std::cos(theta), std::sin(theta), 0, 0;
    partial_u << 0, 0, 1, 0, 0, 0, 0, 1;
    return interval * (x_dot - partial_x.transpose() * x_refer -
                       partial_u.transpose() * u_refer);
  };

  auto steer_rate_cons_matrix_caster =
      [](const TrackerParam &param,
         const TrajectoryTracker::TrajectoryXD &x_refer,
         const TrajectoryTracker::TrajectoryXD &u_refer)
      -> TrajectoryTracker::DMatrix {
    TrajectoryTracker::DMatrix P;
    int wheel_num = 2, state_size = param.state_size_,
        input_size = param.input_size_, horizon = param.horizon_,
        interval = param.interval_, track_width = param.track_width_,
        dist_front_to_rear = param.dist_front_to_rear_,
        qp_state_size = state_size * (horizon + 1) + input_size * horizon;
    P.resize(2 * (horizon - 1), qp_state_size);
    for (size_t i = 0; i < horizon - 1; ++i) {
      TrajectoryTracker::DVector left_cons(qp_state_size);
      left_cons.setZero();
      double v_0 = x_refer.at(i)(3), v_1 = x_refer.at(i + 1)(3),
             omega_0 = u_refer.at(i)(0), omega_1 = u_refer.at(i + 1)(0);
      double g_0 = (2 * dist_front_to_rear * omega_0) /
                   (2 * v_0 - track_width * omega_0),
             g_1 = (2 * dist_front_to_rear * omega_1) /
                   (2 * v_1 - track_width * omega_1);
      double partial_g_v_0 = -(4 * dist_front_to_rear * omega_0) /
                             std::pow(2 * v_0 - track_width * omega_0, 2),
             partial_g_v_1 = -(4 * dist_front_to_rear * omega_1) /
                             std::pow(2 * v_1 - track_width * omega_1, 2),
             partial_g_omege_0 = (4 * dist_front_to_rear * v_0) /
                                 std::pow(2 * v_0 - track_width * omega_0, 2),
             partial_g_omege_1 = (4 * dist_front_to_rear * v_1) /
                                 std::pow(2 * v_1 - track_width * omega_1, 2);
      double partial_h_v_0 =
                 (std::pow(g_1, 2) * partial_g_v_0 + partial_g_v_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_v_1 =
                 (std::pow(g_0, 2) * partial_g_v_1 + partial_g_v_1) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_0 =
                 (std::pow(g_1, 2) * partial_g_omege_0 + partial_g_omege_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_1 =
                 (std::pow(g_0, 2) * partial_g_omege_1 + partial_g_omege_1) /
                 std::pow(1 + g_1 * g_0, 2);
      left_cons(i * state_size + 3) = partial_h_v_0;
      left_cons((i + 1) * state_size + 3) = partial_h_v_1;
      left_cons((horizon + 1) * state_size + i * input_size + 0) =
          partial_g_omege_0;
      left_cons((horizon + 1) * state_size + (i + 1) * input_size + 0) =
          partial_g_omege_1;
      P.row(2 * i) = left_cons.transpose();
    }
    for (size_t i = 0; i < horizon - 1; ++i) {
      TrajectoryTracker::DVector right_cons(qp_state_size);
      right_cons.setZero();
      double v_0 = x_refer.at(i)(3), v_1 = x_refer.at(i + 1)(3),
             omega_0 = u_refer.at(i)(0), omega_1 = u_refer.at(i + 1)(0);
      double g_0 = (2 * dist_front_to_rear * omega_0) /
                   (2 * v_0 + track_width * omega_0),
             g_1 = (2 * dist_front_to_rear * omega_1) /
                   (2 * v_1 + track_width * omega_1);
      double partial_g_v_0 = -(4 * dist_front_to_rear * omega_0) /
                             std::pow(2 * v_0 + track_width * omega_0, 2),
             partial_g_v_1 = -(4 * dist_front_to_rear * omega_1) /
                             std::pow(2 * v_1 + track_width * omega_1, 2),
             partial_g_omege_0 = (4 * dist_front_to_rear * v_0) /
                                 std::pow(2 * v_0 + track_width * omega_0, 2),
             partial_g_omege_1 = (4 * dist_front_to_rear * v_1) /
                                 std::pow(2 * v_1 + track_width * omega_1, 2);
      double partial_h_v_0 =
                 (std::pow(g_1, 2) * partial_g_v_0 + partial_g_v_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_v_1 =
                 (std::pow(g_0, 2) * partial_g_v_1 + partial_g_v_1) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_0 =
                 (std::pow(g_1, 2) * partial_g_omege_0 + partial_g_omege_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_1 =
                 (std::pow(g_0, 2) * partial_g_omege_1 + partial_g_omege_1) /
                 std::pow(1 + g_1 * g_0, 2);
      right_cons(i * state_size + 3) = partial_h_v_0;
      right_cons((i + 1) * state_size + 3) = partial_h_v_1;
      right_cons((horizon + 1) * state_size + i * input_size + 0) =
          partial_g_omege_0;
      right_cons((horizon + 1) * state_size + (i + 1) * input_size + 0) =
          partial_g_omege_1;
      P.row(2 * i + 1) = right_cons.transpose();
    }
    return P;
  };

  auto steer_rate_cons_bound_caster =
      [](const TrackerParam &param,
         const TrajectoryTracker::TrajectoryXD &x_refer,
         const TrajectoryTracker::TrajectoryXD &u_refer)
      -> TrajectoryTracker::DMatrix {
    TrajectoryTracker::DMatrix bound;
    bound.resize(2 * (horizon - 1), 2);
    int horizon = param.horizon_, state_size = param.state_size_,
        input_size = param.input_size_, interval = param.interval_,
        front_wheel_angle_rate_limit = param.front_wheel_angle_rate_limit_,
        qp_state_size = state_size * (horizon + 1) + input_size * horizon;
    for (size_t i = 0; i < horizon - 1; ++i) {
      double v_0 = x_refer.at(i)(3), v_1 = x_refer.at(i + 1)(3),
             omega_0 = u_refer.at(i)(0), omega_1 = u_refer.at(i + 1)(0);
      double g_0 = (2 * dist_front_to_rear * omega_0) /
                   (2 * v_0 - track_width * omega_0),
             g_1 = (2 * dist_front_to_rear * omega_1) /
                   (2 * v_1 - track_width * omega_1);
      double partial_g_v_0 = -(4 * dist_front_to_rear * omega_0) /
                             std::pow(2 * v_0 - track_width * omega_0, 2),
             partial_g_v_1 = -(4 * dist_front_to_rear * omega_1) /
                             std::pow(2 * v_1 - track_width * omega_1, 2),
             partial_g_omege_0 = (4 * dist_front_to_rear * v_0) /
                                 std::pow(2 * v_0 - track_width * omega_0, 2),
             partial_g_omege_1 = (4 * dist_front_to_rear * v_1) /
                                 std::pow(2 * v_1 - track_width * omega_1, 2);
      double partial_h_v_0 =
                 (std::pow(g_1, 2) * partial_g_v_0 + partial_g_v_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_v_1 =
                 (std::pow(g_0, 2) * partial_g_v_1 + partial_g_v_1) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_0 =
                 (std::pow(g_1, 2) * partial_g_omege_0 + partial_g_omege_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_1 =
                 (std::pow(g_0, 2) * partial_g_omege_1 + partial_g_omege_1) /
                 std::pow(1 + g_1 * g_0, 2);
      double k = (g_1 - g_0) / (1 + g_1 * g_0) - partial_h_v_0 * v_0 -
                 partial_h_v_1 * v_1 - partial_h_omege_0 * omega_0 -
                 partial_h_omege_1 * omega_1;
      bound(2 * i, 0) = -std::tan(interval * front_wheel_angle_rate_limit) - k;
      bound(2 * i, 1) = std::tan(interval * front_wheel_angle_rate_limit) - k;
    }
    for (size_t i = 0; i < horizon - 1; ++i) {
      double v_0 = x_refer.at(i)(3), v_1 = x_refer.at(i + 1)(3),
             omega_0 = u_refer.at(i)(0), omega_1 = u_refer.at(i + 1)(0);
      double g_0 = (2 * dist_front_to_rear * omega_0) /
                   (2 * v_0 + track_width * omega_0),
             g_1 = (2 * dist_front_to_rear * omega_1) /
                   (2 * v_1 + track_width * omega_1);
      double partial_g_v_0 = -(4 * dist_front_to_rear * omega_0) /
                             std::pow(2 * v_0 + track_width * omega_0, 2),
             partial_g_v_1 = -(4 * dist_front_to_rear * omega_1) /
                             std::pow(2 * v_1 + track_width * omega_1, 2),
             partial_g_omege_0 = (4 * dist_front_to_rear * v_0) /
                                 std::pow(2 * v_0 + track_width * omega_0, 2),
             partial_g_omege_1 = (4 * dist_front_to_rear * v_1) /
                                 std::pow(2 * v_1 + track_width * omega_1, 2);
      double partial_h_v_0 =
                 (std::pow(g_1, 2) * partial_g_v_0 + partial_g_v_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_v_1 =
                 (std::pow(g_0, 2) * partial_g_v_1 + partial_g_v_1) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_0 =
                 (std::pow(g_1, 2) * partial_g_omege_0 + partial_g_omege_0) /
                 std::pow(1 + g_1 * g_0, 2),
             partial_h_omege_1 =
                 (std::pow(g_0, 2) * partial_g_omege_1 + partial_g_omege_1) /
                 std::pow(1 + g_1 * g_0, 2);
      double k = (g_1 - g_0) / (1 + g_1 * g_0) - partial_h_v_0 * v_0 -
                 partial_h_v_1 * v_1 - partial_h_omege_0 * omega_0 -
                 partial_h_omege_1 * omega_1;
      bound(2 * i + 1, 0) =
          -std::tan(interval * front_wheel_angle_rate_limit) - k;
      bound(2 * i + 1, 1) =
          std::tan(interval * front_wheel_angle_rate_limit) - k;
    }
    return bound;
  };
  // generate reference trajectory
  TrajectoryTracker::Trajectory2D refer_traj;
  refer_traj.reserve(horizon + 1);
  refer_traj.resize(horizon + 1);
  double line_speed = 1.0, radius = 5.0, omega_speed = line_speed / radius;
  for (size_t i = 0; i <= horizon; ++i) {
    double angle = omega_speed * i * interval;
    double x = radius * std::cos(angle), y = radius * std::sin(angle);
    auto &refer_state = refer_traj.at(i);
    refer_state << x, y;
    // refer_state << 5, 5;
  }

  TrajectoryTracker::UniquePtr traj_tracker =
      std::make_unique<TrajectoryTracker>(param);
  traj_tracker->init(Q, R, dynamic_state_matrix_caster,
                     dynamic_input_matrix_caster, dynamic_vector_caster,
                     A_equal, B_equal, K_equal, A_inequal, B_inequal,
                     K_inequal_lb, K_inequal_ub, x_lb, x_ub, u_lb, u_ub,
                     init_state, refer_traj);
  traj_tracker->addUserCustomizedConstraints(steer_rate_cons_matrix_caster,
                                             steer_rate_cons_bound_caster);
  TrajectoryTracker::DVector solution;
  traj_tracker->solve(solution);

  Eigen::IOFormat CleanFmt(4, Eigen::DontAlignCols, ", ", "", "(", ")");
  for (size_t i = 0; i <= horizon; ++i) {
    Eigen::Vector2d x = solution.segment(i * (state_size), 2);
    double dist = (refer_traj.at(i) - x).norm();
    std::cout << std::fixed << std::setprecision(2) << std::setw(4)
              << std::showpos << "time stamp = " << i * interval
              << ", reference state = "
              << refer_traj.at(i).transpose().format(CleanFmt)
              << ", planning state = " << x.transpose().format(CleanFmt)
              << ", error = " << dist << std::endl;
  }
  // traj_tracker->printRefereceStateSeq();
  // traj_tracker->printRefereceInputSeq();
  // traj_tracker->printOsqpMatrices();
}