#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <cmath>

#include "geometry_msgs/Twist.h"
#include "ros/node_handle.h"
#include "ros/timer.h"
#include "tracking_server.h"
#include "trajectory_tracker.h"

using namespace willand_ackermann;

int main(int argc, char** argv) {
  ros::init(argc, argv, "traj_tracking");
  ros::NodeHandle nh;
  int horizon, state_size, input_size;
  double interval, speed_limit, acc_limit, front_wheel_angle_limit,
      front_wheel_angle_rate_limit, track_width, dist_front_to_rear;
  double weight_x_error, weight_y_error, weight_theta_error, weight_v,
      weight_omega;

  nh.param("/traj_tracking/horizon", horizon, 20);
  nh.param("/traj_tracking/interval", interval, 0.05);
  nh.param("/traj_tracking/state_size", state_size, 3);
  nh.param("/traj_tracking/input_size", input_size, 2);
  nh.param("/traj_tracking/speed_limit", speed_limit, 1.5);
  nh.param("/traj_tracking/acc_limit", acc_limit, 1.5);
  nh.param("/traj_tracking/front_wheel_angle_limit", front_wheel_angle_limit,
           M_PI / 2);
  nh.param("/traj_tracking/front_wheel_angle_rate_limit",
           front_wheel_angle_rate_limit, M_PI / 2);
  nh.param("/traj_tracking/track_width", track_width, 0.4);
  nh.param("/traj_tracking/dist_front_to_rear", dist_front_to_rear, 0.4);
  nh.param("/traj_tracking/weight_x_error", weight_x_error, 3333.3);
  nh.param("/traj_tracking/weight_y_error", weight_y_error, 3333.3);
  nh.param("/traj_tracking/weight_theta_error", weight_theta_error, 333.3);
  nh.param("/traj_tracking/weight_v", weight_v, 23.3);
  nh.param("/traj_tracking/weight_omega", weight_omega, 23.3);
  ROS_INFO(
      "horizon: %d, interval: %lf, state_size: %d, input_size: %d, "
      "speed_limit: %lf, acc_limit: %lf, front_wheel_angle_limit: %lf, "
      "front_wheel_angle_rate_limit: %lf, weight_x_error: %lf, weight_y_error: "
      "%lf, weight_theta_error: %lf, weight_v: %lf, weight_omega: %lf, "
      "track_width: %lf, dist_front_to_rear: %lf",
      horizon, interval, state_size, input_size, speed_limit, acc_limit,
      front_wheel_angle_limit, front_wheel_angle_rate_limit, weight_x_error,
      weight_y_error, weight_theta_error, weight_v, weight_omega, track_width,
      dist_front_to_rear);
  TrackerParam param(horizon, interval, state_size, input_size, speed_limit,
                     acc_limit, front_wheel_angle_limit,
                     front_wheel_angle_rate_limit, weight_x_error,
                     weight_y_error, weight_theta_error, weight_v, weight_omega,
                     track_width, dist_front_to_rear);

  constexpr double scale_coef = 0.8;
  double interval_between_points = 1.2 * interval;
  PathGenerator path_generator(interval_between_points);
  TrackingServer tracking_server(param, path_generator);
  tracking_server.init(nh);

  ros::AsyncSpinner spinner(2);  // Use 2 threads
  spinner.start();
  ros::waitForShutdown();

  return 0;
}
