#ifndef SPLINE_OPT_DIS_01
#define SPLINE_OPT_DIS_01

#include "rviz_dis.hpp"
#include "traj_optim/spline_opt.hpp"

struct splineOptDis{
  cubicSplineOpt& cso;
  rviz1DisSimp& r1ds;
  splineOptDis(rviz1DisSimp&_r1ds, cubicSplineOpt& _cso) : 
                 cso(_cso), r1ds(_r1ds){}

  int add_lbfgs_path_points(){
    return r1ds.add_scattered2d(r1ds.id0++, cso.points1, 0.3, 0.0, {0.0,1.0,1.0,1.0});
  }

  int add_lbfgs_path(){
    return r1ds.add_path_strip2d(r1ds.id0++, cso.points1, 0.1, 0.1, {0.0,1.0,0.0,1.0});
  }
};

#endif // SPLINE_OPT_DIS_01
