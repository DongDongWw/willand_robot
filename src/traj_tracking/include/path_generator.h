#include <Eigen/Dense>
#include <boost/signals2/deconstruct.hpp>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

#include "geometry_msgs/Pose2D.h"
#include "ros/ros.h"

namespace simple_ackermann {

class PathGenerator {
 public:
  typedef Eigen::Vector2d Point2D;

 private:
  Point2D start_;
  Point2D end_;
  double interval_;
  std::vector<Point2D> points_;
  Eigen::Matrix<double, 3, 2> coefficients_;

 public:
  PathGenerator() : interval_(0.05){};
  PathGenerator(double interval) : interval_(interval) {}
  std::vector<Point2D> getGlobalPath(const Point2D &start, const Point2D &end,
                                     double theta) {
    start_ = start;
    end_ = end;
    generateCoefficients(theta);
    points_ = std::vector<Point2D>();
    double t = 0.0;
    constexpr int max_points_num = 1000;
    while (std::abs(t - 1.0) > 1e-6 && points_.size() < max_points_num) {
      points_.push_back(getPoint(t));
      double next_t = findParameterForArcLength(interval_, t);
      t = next_t;
    }
    points_.push_back(getPoint(t));
    return points_;
  }
  std::vector<Point2D> getGlobalPath(double radius) {
    points_.clear();
    double radius_abs = std::abs(radius);
    double delta_theta = interval_ / radius_abs;
    // positive radius means turning left, anti-clockwise
    for (double theta = 0; theta < 2 * M_PI; theta += delta_theta) {
      double x = radius_abs * std::cos(theta - M_PI / 2);
      double y = radius_abs * std::sin(theta - M_PI / 2) + radius_abs;
      if (radius < 0) {
        y = -y;
      }
      points_.push_back(Point2D(x, y));
    }
    return points_;
  }
  // TODO: add dubins path
  // std::vector<Point2D> getGlobalPath(const Point2D &target_point)
  std::vector<Point2D> generateReferenceTrajectory(const Point2D &p,
                                                   int number_of_points) {
    if (points_.empty()) {
      return std::vector<Point2D>();
    }
    double min_dist = std::numeric_limits<double>::max();
    int min_dist_idx = 0;
    for (size_t i = 0; i < points_.size(); ++i) {
      double dist = (points_.at(i) - p).norm();
      if (dist < min_dist) {
        min_dist = dist;
        min_dist_idx = i;
      }
    }
    std::vector<Point2D> ref_traj(number_of_points, points_.back());
    for (size_t i = min_dist_idx;
         i < min_dist_idx + number_of_points && i < points_.size(); ++i) {
      ref_traj.at(i - min_dist_idx) = points_.at(i);
    }
    return ref_traj;
  }
  double getDistOffset(const Point2D &p) {
    if (points_.empty()) {
      return 0.0;
    }
    double min_dist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < points_.size(); ++i) {
      double dist = (points_.at(i) - p).norm();
      if (dist < min_dist) {
        min_dist = dist;
      }
    }
    return min_dist;
  }

 private:
  void generateCoefficients(double theta) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(2.0, 4.0);

    coefficients_(0, 0) = start_(0);
    coefficients_(0, 1) = start_(1);
    if (std::abs(theta - M_PI / 2) < 1e-6) {
      coefficients_(1, 0) = 0;
      coefficients_(1, 1) = dis(gen);
    } else if (std::abs(theta + M_PI / 2) < 1e-6) {
      coefficients_(1, 0) = 0;
      coefficients_(1, 1) = -dis(gen);
    } else {
      double slope = std::tan(theta);
      if (theta > -M_PI / 2 && theta < M_PI / 2) {
        coefficients_(1, 0) = dis(gen);
        coefficients_(1, 1) = slope * coefficients_(1, 0);
      } else {
        coefficients_(1, 0) = -dis(gen);
        coefficients_(1, 1) = slope * coefficients_(1, 0);
      }
    }
    coefficients_(2, 0) = end_(0) - coefficients_(0, 0) - coefficients_(1, 0);
    coefficients_(2, 1) = end_(1) - coefficients_(0, 1) - coefficients_(1, 1);
  }

  Point2D getPoint(double t) const {
    Eigen::VectorXd powers(3);
    powers << 1, t, t * t;
    return coefficients_.transpose() * powers;
  }

  Point2D getDerivative(double t) const {
    Eigen::VectorXd powers(3);
    powers << 0, 1, 2 * t;
    return coefficients_.transpose() * powers;
  }

  double calculateArcLength(double t0, double t1, int num_samples = 10) const {
    double length = 0.0;
    double dt = (t1 - t0) / num_samples;
    for (int i = 0; i < num_samples; ++i) {
      double t = t0 + i * dt;
      Point2D derivative = getDerivative(t);
      length += derivative.norm() * dt;
    }
    return length;
  }

  double findParameterForArcLength(double target_length, double start,
                                   int num_samples = 10) const {
    double low = start, high = 1.0;
    double length = calculateArcLength(start, 1.0, num_samples);
    if (length < target_length) {
      return 1.0;
    }
    while (high - low > 1e-6) {
      double mid = (low + high) / 2.0;
      double length = calculateArcLength(start, mid, num_samples);
      if (length < target_length) {
        low = mid;
      } else {
        high = mid;
      }
    }
    return (low + high) / 2.0;
  }
};
};  // namespace simple_ackermann