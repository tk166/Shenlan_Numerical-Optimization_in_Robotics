// class to draw in rviz 1
// V0.0.1 20220825 tkalpha
#ifndef RVIZ_DIS_01
#define RVIZ_DIS_01
#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>

#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <string>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;

struct colorSimp {
  double r;
  double g;
  double b;
  double a;
  colorSimp(double r0 = 0.5, double g0 = 0.6, double b0 = 0.5, double a0 = 1.0) : r(r0), g(g0), b(b0), a(a0) {}
  ~colorSimp() {}
};

struct rviz1DisSimp {
  ros::NodeHandle n;
  ros::Publisher pub;
  ros::Rate r;
  std::string p_name;
  std::string f_name;
  std::string ns;
  visualization_msgs::MarkerArray ma;
  int id0;
  rviz1DisSimp(ros::NodeHandle& nh, const std::string publisher_name, const std::string frame_name, const std::string name_space)
      : n(nh), r(30), p_name(publisher_name), f_name(frame_name), ns(name_space), id0(0) {
    pub = n.advertise<visualization_msgs::MarkerArray>(publisher_name, 10);
    while (pub.getNumSubscribers() == 0) {
      r.sleep();
    }
  }
  virtual ~rviz1DisSimp() {}

  visualization_msgs::Marker get_default_marker(int id) {
    visualization_msgs::Marker m;
    m.header.frame_id = f_name;
    m.header.stamp = ros::Time::now();
    m.ns = ns;
    m.action = visualization_msgs::Marker::ADD;
    m.id = id;
    m.type = visualization_msgs::Marker::LINE_STRIP;
    m.pose.orientation.w = 1.0;
    m.scale.x = 0.2;
    m.color.r = m.color.g = m.color.b = m.color.a = 1.0;
    return m;
  }

  int add_convex_poly_flat(int id, Mat& vertex, double height, colorSimp c = colorSimp()) {
    visualization_msgs::Marker m = get_default_marker(id);
    m.type = visualization_msgs::Marker::TRIANGLE_LIST;
    m.scale.x = m.scale.y = m.scale.z = 1.0;
    geometry_msgs::Point pl;
    std_msgs::ColorRGBA c1;
    c1.r = c.r;
    c1.g = c.g;
    c1.b = c.b;
    c1.a = c.a;
    m.points.clear();
    m.colors.clear();
    int num = vertex.rows();
    // add the display of top and bottom surfaces
    for (int k = 1; k <= num - 2; ++k) {
      pl.z = 0.0;
      pl.x = vertex(0, 0);
      pl.y = vertex(0, 1);
      m.points.push_back(pl);
      pl.x = vertex(k, 0);
      pl.y = vertex(k, 1);
      m.points.push_back(pl);
      pl.x = vertex(k + 1, 0);
      pl.y = vertex(k + 1, 1);
      m.points.push_back(pl);
      m.colors.push_back(c1);
      pl.z = height;
      pl.x = vertex(0, 0);
      pl.y = vertex(0, 1);
      m.points.push_back(pl);
      pl.x = vertex(k, 0);
      pl.y = vertex(k, 1);
      m.points.push_back(pl);
      pl.x = vertex(k + 1, 0);
      pl.y = vertex(k + 1, 1);
      m.points.push_back(pl);
      m.colors.push_back(c1);
    }
    std_msgs::ColorRGBA c2;
    c2.r = c.r * 0.8;
    c2.g = c.g * 0.8;
    c2.b = c.b * 0.8;
    c2.a = c.a;
    // add the display of surfaces on the side
    int j = num - 1;
    for (int k = 0; k <= num - 1; ++k) {
      pl.z = height;
      pl.x = vertex(k, 0);
      pl.y = vertex(k, 1);
      m.points.push_back(pl);
      pl.z = 0.0;
      m.points.push_back(pl);
      pl.x = vertex(j, 0);
      pl.y = vertex(j, 1);
      m.points.push_back(pl);
      m.colors.push_back(c2);
      m.points.push_back(pl);
      pl.z = height;
      m.points.push_back(pl);
      pl.x = vertex(k, 0);
      pl.y = vertex(k, 1);
      m.points.push_back(pl);
      m.colors.push_back(c2);
      j = k;
    }
    ma.markers.push_back(m);
    return 0;
  }

  int add_cylinder(int id, double x0, double y0, double r0, double height, colorSimp c = colorSimp()) {
    visualization_msgs::Marker m = get_default_marker(id);
    m.type = visualization_msgs::Marker::CYLINDER;
    m.color.r = c.r;
    m.color.g = c.g;
    m.color.b = c.b;
    m.color.a = c.a;
    m.pose.position.x = x0;
    m.pose.position.y = y0;
    m.scale.x = m.scale.y = r0 * 2;
    m.scale.z = height;
    ma.markers.push_back(m);
    return 0;
  }

  int add_path_strip2d(int id, Mat& path, double width, double height, colorSimp c = colorSimp()) {
    visualization_msgs::Marker m = get_default_marker(id);
    m.type = visualization_msgs::Marker::LINE_STRIP;
    m.color.r = c.r;
    m.color.g = c.g;
    m.color.b = c.b;
    m.color.a = c.a;
    m.points.clear();
    m.scale.x = width;
    geometry_msgs::Point pl;
    pl.z = height;
    for (int k = 0; k < path.rows(); ++k) {
      pl.x = path(k, 0);
      pl.y = path(k, 1);
      m.points.push_back(pl);
    }
    ma.markers.push_back(m);
    return 0;
  }

  int add_scattered2d(int id, Mat& points, double diameter, double height, colorSimp c = colorSimp()) {
    visualization_msgs::Marker m = get_default_marker(id);
    m.type = visualization_msgs::Marker::SPHERE_LIST;
    m.color.r = c.r;
    m.color.g = c.g;
    m.color.b = c.b;
    m.color.a = c.a;
    m.points.clear();
    m.scale.x = m.scale.y = m.scale.z = diameter;
    geometry_msgs::Point pl;
    pl.z = height;
    for (int k = 0; k < points.rows(); ++k) {
      pl.x = points(k, 0);
      pl.y = points(k, 1);
      m.points.push_back(pl);
    }
    ma.markers.push_back(m);
    return 0;
  }

  int add_colored_path_strip2d(int id, Mat& path, Vec& val, double val_min, double val_max, double width,
                               double height, colorSimp c_min = colorSimp(0.8, 0.2, 0.2, 1.0), colorSimp c_max = colorSimp(1.0, 1.0, 0.2, 1.0)) {
    visualization_msgs::Marker m = get_default_marker(id);
    m.type = visualization_msgs::Marker::LINE_STRIP;
    m.color.r = 1.0;
    m.color.g = 1.0;
    m.color.b = 1.0;
    m.color.a = 1.0;
    m.points.clear();
    m.scale.x = width;
    geometry_msgs::Point pl;
    std_msgs::ColorRGBA cl;
    pl.z = height;
    for (int k = 0; k < path.rows() - 1; ++k) {
      pl.x = path(k, 0);
      pl.y = path(k, 1);
      double cl_pos = std::min(1.0, std::max(0.0, (val(k) - val_min) / (val_max - val_min)));
      cl.r = c_min.r + (c_max.r - c_min.r) * cl_pos;
      cl.g = c_min.g + (c_max.g - c_min.g) * cl_pos;
      cl.b = c_min.b + (c_max.b - c_min.b) * cl_pos;
      cl.a = 1.0;
      m.points.push_back(pl);
      m.colors.push_back(cl);
    }
    ma.markers.push_back(m);
    return 0;
  }

  int add_arrow(int id, double x0, double y0, double z0, double x1, double y1, double z1, double size, colorSimp c = colorSimp(1.0, 0, 0, 1.0)) {
    visualization_msgs::Marker m = get_default_marker(id);
    m.type = visualization_msgs::Marker::ARROW;
    m.color.r = c.r;
    m.color.g = c.g;
    m.color.b = c.b;
    m.color.a = c.a;
    m.scale.x = size * 0.5;
    m.scale.y = size;
    m.scale.z = size;
    geometry_msgs::Point pl;
    pl.x = x0;
    pl.y = y0;
    pl.z = z0;
    m.points.push_back(pl);
    pl.x = x1;
    pl.y = y1;
    pl.z = z1;
    m.points.push_back(pl);
    ma.markers.push_back(m);
    return 0;
  }

  int clear() {
    ma.markers.clear();
    id0 = 0;
    return 0;
  }
  int send() {
    if (ros::ok()) {
      pub.publish(ma);
      pub.publish(ma);
    } else {
      std::cout << "ros NOT ok." << std::endl;
    }
    return 0;
  }
};

#endif  // RVIZ_DIS_01