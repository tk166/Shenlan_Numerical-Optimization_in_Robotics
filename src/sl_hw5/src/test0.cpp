#include <iostream>

#include "ros1_disp/rviz_dis.hpp"
#include "scene/pavel_valtr.hpp"
#include "sdqp/sdqp.hpp"

int main(int argc, char** argv) {
  ros::init(argc, argv, "sl_hw5");
  ros::NodeHandle nh;
  // initialize a polygon
  rviz1DisSimp rd0(nh, "sl_hw5_array", "map", "sl_hw5_test");
  Mat poly1(4, 2);
  poly1 << 0, 0, 1, 0, 2, 1, 0, 1;
  randPoly rp;
  rp.generate(8, 5.6, poly1);

  // formulate the constraints
  double x0 = -5, y0 = -4;
  int num = poly1.rows();
  Eigen::Matrix<double, 2, 2> Q(Mat::Identity(2, 2));
  Eigen::Matrix<double, 2, 1> c;
  Eigen::Matrix<double, 2, 1> x;
  c << -2 * x0, -2 * y0;
  Mat A(num, 2);
  Vec b(num);
  int j = num - 1;
  for (int k = 0; k < num; ++k) {
    A(k, 0) = -(poly1(j, 1) - poly1(k, 1));
    A(k, 1) = -(poly1(k, 0) - poly1(j, 0));
    b(k) = (A(k, 0) * poly1(j, 0) + A(k, 1) * poly1(j, 1));
    j = k;
  }
  double minobj = sdqp::sdqp<2>(Q, c, A, b, x);
  std::cout << Q << "\n\n"
            << c << "\n\n"
            << A << "\n\n"
            << b << std::endl;
  std::cout << x.transpose() << "    " << minobj << std::endl;

  // draw out the scene
  rd0.add_convex_poly_flat(0, poly1, 0.1);
  rd0.add_arrow(1, x0, y0, 0, x(0), x(1), 0, 0.2);
  rd0.send();
}

// #include <ros/ros.h>
// #include <visualization_msgs/Marker.h>

// #include <cmath>

// int main(int argc, char** argv) {
//   ros::init(argc, argv, "points_and_lines");
//   ros::NodeHandle n;
//   ros::Publisher marker_pub =
//       n.advertise<visualization_msgs::Marker>("sl_hw5_array", 10);

//   ros::Rate r(30);

//   float f = 0.0;
//   while (ros::ok()) {
//     // %Tag(MARKER_INIT)%
//     visualization_msgs::Marker points, line_strip, line_list;
//     points.header.frame_id = line_strip.header.frame_id =
//         line_list.header.frame_id = "map";
//     points.header.stamp = line_strip.header.stamp = line_list.header.stamp =
//         ros::Time::now();
//     points.ns = line_strip.ns = line_list.ns = "points_and_lines";
//     points.action = line_strip.action = line_list.action =
//         visualization_msgs::Marker::ADD;
//     points.pose.orientation.w = line_strip.pose.orientation.w =
//         line_list.pose.orientation.w = 1.0;
//     // %EndTag(MARKER_INIT)%

//     // %Tag(ID)%
//     points.id = 0;
//     line_strip.id = 1;
//     line_list.id = 2;
//     // %EndTag(ID)%

//     // %Tag(TYPE)%
//     points.type = visualization_msgs::Marker::POINTS;
//     line_strip.type = visualization_msgs::Marker::LINE_STRIP;
//     line_list.type = visualization_msgs::Marker::LINE_LIST;
//     // %EndTag(TYPE)%

//     // %Tag(SCALE)%
//     // POINTS markers use x and y scale for width/height respectively
//     points.scale.x = 0.2;
//     points.scale.y = 0.2;

//     // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the
//     // line width
//     line_strip.scale.x = 0.1;
//     line_list.scale.x = 0.1;
//     // %EndTag(SCALE)%

//     // %Tag(COLOR)%
//     // Points are green
//     points.color.g = 1.0f;
//     points.color.a = 1.0;

//     // Line strip is blue
//     line_strip.color.b = 1.0;
//     line_strip.color.a = 1.0;

//     // Line list is red
//     line_list.color.r = 1.0;
//     line_list.color.a = 1.0;
//     // %EndTag(COLOR)%

//     // %Tag(HELIX)%
//     // Create the vertices for the points and lines
//     for (uint32_t i = 0; i < 100; ++i) {
//       float y = 5 * sin(f + i / 100.0f * 2 * M_PI);
//       float z = 5 * cos(f + i / 100.0f * 2 * M_PI);

//       geometry_msgs::Point p;
//       p.x = (int32_t)i - 50;
//       p.y = y;
//       p.z = z;

//       points.points.push_back(p);
//       line_strip.points.push_back(p);

//       // The line list needs two points for each line
//       line_list.points.push_back(p);
//       p.z += 1.0;
//       line_list.points.push_back(p);
//     }
//     // %EndTag(HELIX)%

//     marker_pub.publish(points);
//     marker_pub.publish(line_strip);
//     marker_pub.publish(line_list);

//     r.sleep();

//     f += 0.04;
//   }
// }
// // %EndTag(FULLTEXT)%
