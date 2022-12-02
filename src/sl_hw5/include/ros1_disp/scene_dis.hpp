#ifndef SCENE_DISP_01
#define SCENE_DISP_01

#include "rviz_dis.hpp"
#include "scene/scene_astar.hpp"

struct sceneAstarDisp {
  sceneAstar& sa;
  rviz1DisSimp& r1ds;
  sceneAstarDisp(rviz1DisSimp&_r1ds, sceneAstar& _sa) : 
                 sa(_sa), r1ds(_r1ds){}
  int clear(){
    return r1ds.clear();
  }
  int send(){
    return r1ds.send();
  }
  int add_scene_grid(){
    Mat points(sa.num_x*sa.num_y, 2);
    for(size_t i = 0; i<sa.num_x; ++i){
      for(size_t j = 0; j<sa.num_y; ++j){
        points(i+sa.num_x*j,0) = sa.map0[i][j].pos(0);
        points(i+sa.num_x*j,1) = sa.map0[i][j].pos(1);
      }
    }
    return r1ds.add_scattered2d(r1ds.id0++, points, 0.1, 0.0, {1.0,1.0,1.0,1.0});
  }

  int add_pos_goal(){
    Mat points(2, 2);
    points(0, 0) = sa.pos_start(0);
    points(0, 1) = sa.pos_start(1);
    points(1, 0) = sa.pos_goal(0);
    points(1, 1) = sa.pos_goal(1);
    return r1ds.add_scattered2d(r1ds.id0++, points, 1.0, 0.0, {0.0,1.0,0.0,1.0});
  }

  int add_obstacles(){
    for(auto& ob_ptr : sa.obs){
      if(ob_ptr->type == 0){
        auto ob_ptr1 = (cylinderObstacle*) ob_ptr;
        r1ds.add_cylinder(r1ds.id0++, ob_ptr1->x, ob_ptr1->y, ob_ptr1->r, 0.2);
      } else {
        auto ob_ptr1 = (polyObstacle*) ob_ptr;
        r1ds.add_convex_poly_flat(r1ds.id0++, ob_ptr1->poly, 0.2);
      }
    }
    return 0;
  }

  int add_astar_path(){
    return r1ds.add_path_strip2d(r1ds.id0++, sa.astar_path, 0.1, 0.1, {0.0,0.7,0.0,1.0});
  }
};

#endif  // SCENE_DISP_01