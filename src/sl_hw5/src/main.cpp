#include <iostream>

#include "ros1_disp/scene_dis.hpp"
#include "ros1_disp/spline_opt_dis.hpp"
#include "ros1_disp/topp_dis.hpp"
#include "scene/scene_astar.hpp"
#include "traj_optim/spline_opt.hpp"
#include "traj_optim/conic_alm_topp.hpp"

int main(int argc, char** argv) {
  ros::init(argc, argv, "sl_hw5");
  ros::NodeHandle nh;
  // initialize a scene
  sceneAstar sa(-50., -25., 0., 100., 50., 0.5, 3.6, 12.9, 3, 12);
  std::cout << "[[[ scene generation OK!]]]" << std::endl;
  sa.random_obstacles(26);
  std::cout << "[[[ obstacles generation OK!]]]" << std::endl;
  // std:: cout <<"{{{"<<sa.is_free(0.0, 0.0) <<std::endl;;
  int search_state = sa.random_astar_path();
  std::cout << "[[[ Astart search state: " << search_state << " ]]]" << std::endl;

  // path optimization
  cubicSplineOpt cso(sa, sa.astar_path, 0.9, 5., -18.);
  std::cout << "[[[ initial cubic spline OK!]]]" << std::endl;
  cso.path_opt_lbfgs();
  std::cout << "[[[ L-BFGS path optimization OK!]]]" << std::endl;

  // trajectory time optimization
  cso.topp_prepare(2);
  conicALMTOPP2 topp2(cso.s1, cso.q1, cso.qv1, cso.qa1, 1.0, 10.0, 0.0, 0.0);
  Vec a,b,c,d;
  double result;
  int topp_ret = topp2.solve(result, a, b, c, d);
  std::cout << "\n===============\n[Solving status:] " << topp_ret << "\n";
  std::cout << "[Minimum value:] " << result << "\n";
  std::cout << "[a:] " << a.transpose() << std::endl;
  std::cout << "[b:] " << b.transpose() << std::endl;
  std::cout << "[c:] " << c.transpose() << std::endl;
  std::cout << "[d:] " << d.transpose() << std::endl;


  
  // draw out the scene
  rviz1DisSimp disp_rviz(nh, "sl_hw5_array", "map", "sl_hw5_test");
  sceneAstarDisp disp_sa(disp_rviz, sa);
  std::cout << "[[[ disp_sa initialization OK!]]]" << std::endl;
  disp_sa.add_scene_grid();
  disp_sa.add_pos_goal();
  std::cout << "[[[ disp_sa add grid OK!]]]" << std::endl;
  disp_sa.add_obstacles();
  std::cout << "[[[ disp_sa add obstacles OK!]]]" << std::endl;
  disp_sa.add_astar_path();
  std::cout << "[[[ disp_sa add Astar path OK!]]]" << std::endl;

  // splineOptDis disp_cso(disp_rviz, cso);
  // disp_cso.add_lbfgs_path_points();
  // disp_cso.add_lbfgs_path();
  // std::cout << "[[[ disp_sa add L-BFGS path OK!]]]" << std::endl;

  conicAlmToppDis disp_topp(disp_rviz, topp2);
  disp_topp.add_topp_trajectory_points(0.8, 6.8);
  std::cout << "[[[ disp_sa add TOPP trajectory OK!]]]" << std::endl;
  disp_rviz.send();
  disp_rviz.send();
}
