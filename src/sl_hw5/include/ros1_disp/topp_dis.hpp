#ifndef TOPP_DIS_01
#define TOPP_DIS_01

#include <fstream>

#include "rviz_dis.hpp"
#include "traj_optim/conic_alm_topp.hpp"

struct conicAlmToppDis {
  conicALMTOPP2& topp2;
  rviz1DisSimp& r1ds;
  conicAlmToppDis(rviz1DisSimp& _r1ds, conicALMTOPP2& _topp2) : topp2(_topp2), r1ds(_r1ds) {}

  int add_topp_trajectory_points(double v_min, double v_max) {
    Vec a, b, c, d;
    topp2.get_result(a, b, c, d);
    Vec v = b.cwiseSqrt();
    std::ofstream file;
    file.open("/home/tku/1/dtt1_q.txt");
    for(int k = 0; k < (int)topp2.q.rows(); ++k){
      file << topp2.q(k,0)<<", "<<topp2.q(k,1)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_qv.txt");
    for(int k = 0; k < (int)topp2.qv.rows(); ++k){
      file << topp2.qv(k,0)<<", "<<topp2.qv(k,1)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_qa.txt");
    for(int k = 0; k < (int)topp2.qa.rows(); ++k){
      file << topp2.qa(k,0)<<", "<<topp2.qa(k,1)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_s.txt");
    for(int k = 0; k < (int)topp2.s.size(); ++k){
      file << topp2.s(k)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_a.txt");
    for(int k = 0; k < (int)a.rows(); ++k){
      for(int j = 0; j < (int)a.cols()-1; ++j){
        file << a(k,j)<<", ";
      }
      file << a(k,a.cols()-1)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_b.txt");
    for(int k = 0; k < (int)b.rows(); ++k){
      for(int j = 0; j < (int)b.cols()-1; ++j){
        file << b(k,j)<<", ";
      }
      file << b(k,b.cols()-1)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_c.txt");
    for(int k = 0; k < (int)c.rows(); ++k){
      for(int j = 0; j < (int)c.cols()-1; ++j){
        file << c(k,j)<<", ";
      }
      file << c(k,c.cols()-1)<<"\n";
    }
    file.close();
    file.open("/home/tku/1/dtt1_d.txt");
    for(int k = 0; k < (int)d.rows(); ++k){
      for(int j = 0; j < (int)d.cols()-1; ++j){
        file << d(k,j)<<", ";
      }
      file << d(k,d.cols()-1)<<"\n";
    }
    file.close();
    return r1ds.add_colored_path_strip2d(r1ds.id0++, topp2.q, v, v_min, v_max, 0.25, 0.1);
  }
};

#endif  // TOPP_DIS_01