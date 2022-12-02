#ifndef SCENE_BASE_01
#define SCENE_BASE_01

#include <eigen3/Eigen/Core>
#include <map>
#include <memory>
#include <chrono>

typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;

// Base class of a scene that can provide artifital distance field information
class sceneBase {
 public:
  virtual double dist_field(const double x1, const double y1, Vec* grad = nullptr) = 0;
  virtual ~sceneBase() {}
};

#endif  // SCENE_BASE_01