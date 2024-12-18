#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/Odometry.h>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>

#include <cmath>

#include "geometry_msgs/Twist.h"
#include "mpc_tracker.h"
#include "ros/node_handle.h"
#include "ros/timer.h"
#include "tracking_server.h"

using namespace simple_ackermann;

int main(int argc, char** argv) {
  ros::init(argc, argv, "traj_tracking");
  ros::NodeHandle nh;
  int horizon, state_size, input_size;
  double interval, speed_limit, acc_limit, front_wheel_angle_limit,
      front_wheel_angle_rate_limit, track_width, dist_front_to_rear;
  double weight_x_error, weight_y_error, weight_theta_error, weight_v,
      weight_omega;

  TrackingServer tracking_server(nh);
  //   ros::AsyncSpinner spinner(2);  // Use 2 threads
  //   spinner.start();
  //   ros::waitForShutdown();
  ros::spin();
  return 0;
}
