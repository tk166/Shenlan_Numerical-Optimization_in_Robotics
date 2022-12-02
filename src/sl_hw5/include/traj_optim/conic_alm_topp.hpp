// Reconstruction of the origin conic ALM TOPP algorithm
// with the simplication of unnecessary matrix operations
// and minor fixes in ALM algorithm steps
// V0.1.5 20220831, tkalpha

#ifndef CONIC_ALM_TOPP_01
#define CONIC_ALM_TOPP_01

#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCore>
#include <iostream>
#include <vector>

#include "lbfgs/lbfgs.hpp"

typedef Eigen::Matrix3d Mat3;
typedef Eigen::Vector3d Vec3;
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::SparseMatrix<double> SpMat;

struct conicALMTOPP2 {
  // path information
  Vec& s;
  Mat& q;
  Mat& qv;
  Mat& qa;
  int mode;
  double a_max;
  double v_max;
  double v_start;
  double v_end;

  // variables for conic alm algorithm
  Vec c;
  std::vector<SpMat> As;
  std::vector<Vec3> bs;
  SpMat Aieq;
  SpMat Aieq_transpose;
  Vec bieq;
  SpMat G;
  SpMat G_transpose;
  Vec h;
  double value;
  Vec x0;
  int pts, dim;
  int idx_a, idx_b, idx_c, idx_d;
  int dim_a, dim_b, dim_c, dim_d;

  // ALM algorithm variables
  std::vector<Vec3> mius;
  Vec eta;
  Vec lambda;
  Vec zeros_ieq;

  // L-BFGS parameters
  lbfgs::lbfgs_parameter_t pl2;

  // ALM algorithm parameters
  struct alm_param {
    double alm_rho = 1.0;
    double alm_beta = 1e3;
    double alm_gamma = 1.0;
    double alm_xi = 0.9;
    double alm_xi_min = 1e-8;
    double alm_epsilon_cons = 1e-3;
    double alm_epsilon_prec = 1e-3;
    int alm_iter_max = 80;
    int alm_iter_display = 1;
  } pa2;

  conicALMTOPP2(Vec& _s, Mat& _q, Mat& _qv, Mat& _qa, double _a_max,
                double _v_max, double _v_start, double _v_end)
      : s(_s),
        q(_q),
        qv(_qv),
        qa(_qa),
        mode(0),
        a_max(_a_max),
        v_max(_v_max),
        v_start(_v_start),
        v_end(_v_end),
        value(0),
        pl2(lbfgs::lbfgs_parameter_t()) {
    generate_topp_matrices();
    generate_alm_variables();
  }

  // basic matrices for topp socp problem
  int generate_topp_matrices() {
    pts = s.size();
    dim = 4 * pts - 2;
    idx_a = 3 * pts - 1;
    idx_b = 0;
    idx_c = pts;
    idx_d = 2 * pts;
    dim_a = pts - 1;
    dim_b = pts;
    dim_c = pts;
    dim_d = pts - 1;
    // goal
    x0 = Vec::Zero(dim);
    c = Vec::Zero(dim);
    for (int k = 0; k < dim_d; ++k) {
      c(idx_d + k) = 2 * (s(k + 1) - s(k));
    }
    // cones
    As.clear();
    bs.clear();
    for (int k = 0; k < dim_d; ++k) {
      SpMat A(3, dim);
      A.setZero();
      Vec3 b(Vec3::Zero());
      b(1) = 2;
      A.insert(0, idx_c + k) = A.insert(0, idx_c + k + 1) = 1;
      A.insert(0, idx_d + k) = 1;
      A.insert(2, idx_c + k) = A.insert(2, idx_c + k + 1) = 1;
      A.insert(2, idx_d + k) = -1;
      As.push_back(A);
      bs.push_back(b);
    }
    for (int k = 0; k < dim_b; ++k) {
      SpMat A(3, dim);
      A.setZero();
      Vec3 b(Vec3::Zero());
      b(0) = 1;
      b(2) = -1;
      A.insert(1, idx_c + k) = 2;
      A.insert(0, idx_b + k) = A.insert(2, idx_b + k) = 1;
      As.push_back(A);
      bs.push_back(b);
    }
    // inequality
    int num_ieq = 3 * dim_b + 4 * dim_a;
    Aieq.resize(num_ieq, dim);
    Aieq.setZero();
    bieq = Vec::Zero(num_ieq);
    int pos = 0;
    for (int k = 0; k < dim_b; ++k) {
      Aieq.insert(pos + k, idx_b + k) = -1;
    }
    pos += dim_b;
    for (int k = 0; k < dim_b; ++k) {
      Aieq.insert(pos + k, idx_b + k) = std::pow(qv(k, 0), 2);
      bieq(pos + k) = std::pow(v_max, 2);
    }
    pos += dim_b;
    for (int k = 0; k < dim_b; ++k) {
      Aieq.insert(pos + k, idx_b + k) = std::pow(qv(k, 1), 2);
      bieq(pos + k) = std::pow(v_max, 2);
    }
    pos += dim_b;
    for (int k = 0; k < dim_a; ++k) {
      Aieq.insert(pos + k, idx_b + k) = qa(k, 0);
      Aieq.insert(pos + k, idx_a + k) = qv(k, 0);
      bieq(pos + k) = a_max;
    }
    pos += dim_a;
    for (int k = 0; k < dim_a; ++k) {
      Aieq.insert(pos + k, idx_b + k) = -qa(k, 0);
      Aieq.insert(pos + k, idx_a + k) = -qv(k, 0);
      bieq(pos + k) = a_max;
    }
    pos += dim_a;
    for (int k = 0; k < dim_a; ++k) {
      Aieq.insert(pos + k, idx_b + k) = qa(k, 1);
      Aieq.insert(pos + k, idx_a + k) = qv(k, 1);
      bieq(pos + k) = a_max;
    }
    pos += dim_a;
    for (int k = 0; k < dim_a; ++k) {
      Aieq.insert(pos + k, idx_b + k) = -qa(k, 1);
      Aieq.insert(pos + k, idx_a + k) = -qv(k, 1);
      bieq(pos + k) = a_max;
    }
    // equality
    int num_eq = dim_a + 2;
    G.resize(num_eq, dim);
    G.setZero();
    h = Vec::Zero(num_eq);
    pos = 0;
    for (int k = 0; k < dim_a; ++k) {
      G.insert(pos + k, idx_b + k) = 1;
      G.insert(pos + k, idx_b + k + 1) = -1;
      G.insert(pos + k, idx_a + k) = 2 * (s(k + 1) - s(k));
    }
    pos += dim_a;
    G.insert(pos, idx_b) = qv.row(0).squaredNorm();
    h(pos) = std::pow(v_start, 2);
    G.insert(pos + 1, idx_b + dim_b - 1) = qv.row(dim_b - 1).squaredNorm();
    h(pos + 1) = std::pow(v_end, 2);
    // simplify
    Aieq_transpose = Aieq.transpose();
    G_transpose = G.transpose();
    return 0;
  }

  // initialize ALM variables
  int generate_alm_variables() {
    mius.clear();
    int n_cones = As.size();
    for (int k = 0; k < n_cones; ++k) {
      mius.push_back(Vec3::Zero());
    }
    eta = Vec::Zero(Aieq.rows());
    lambda = Vec::Zero(G.rows());
    zeros_ieq = Vec::Zero(Aieq.rows());
    return 0;
  }

  // function to generate projection on 3d second order cone
  int get_conic3_prj(const Vec3& v, Vec3& prj) {
    double v0 = v(0);
    double v1_norm = std::sqrt(v(1) * v(1) + v(2) * v(2));
    if (v0 <= -v1_norm) {
      prj = Vec3::Zero();
    } else if (v0 >= v1_norm) {
      prj = v;
    } else {
      double c = (v0 + v1_norm) / v1_norm / 2;
      prj(0) = c * v1_norm;
      prj(1) = c * v(1);
      prj(2) = c * v(2);
    }
    return 0;
  }

  // objective function of conic ALM
  static double obj_fcn_alm(void* ptr_obj, const Vec& x, Vec& grad) {
    conicALMTOPP2& p = *(reinterpret_cast<conicALMTOPP2*>(ptr_obj));

    double rho = p.pa2.alm_rho;
    Vec& c = p.c;
    std::vector<SpMat>& As = p.As;
    std::vector<Vec3>& bs = p.bs;
    SpMat& Aieq = p.Aieq;
    SpMat& Aieq_transpose = p.Aieq_transpose;
    Vec& bieq = p.bieq;
    SpMat& G = p.G;
    SpMat& G_transpose = p.G_transpose;
    Vec& h = p.h;

    std::vector<Vec3>& mius = p.mius;
    Vec& eta = p.eta;
    Vec& lambda = p.lambda;
    Vec& zeros_ieq = p.zeros_ieq;

    double value = c.dot(x);
    grad = c;

    Vec3 prj;
    double value_cone = 0;
    int n_cones = As.size();
    for (int k = 0; k < n_cones; ++k) {
      p.get_conic3_prj((1.0 / rho) * mius[k] - As[k] * x - bs[k], prj);
      value_cone += 0.5 * rho * prj.squaredNorm();
      grad -= rho * As[k].transpose() * prj;
    }

    Vec res_ieq = ((1.0 / rho) * eta + Aieq * x - bieq).cwiseMax(zeros_ieq);
    double value_ieq = 0.5 * rho * res_ieq.squaredNorm();
    grad += rho * Aieq_transpose * res_ieq;

    Vec res_eq = (1.0 / rho) * lambda + G * x - h;
    double value_eq = 0.5 * rho * res_eq.squaredNorm();
    grad += rho * G_transpose * res_eq;

    value = value + value_cone + value_ieq + value_eq;
    return value;
  }

  // get ALM error
  double get_alm_error(const Vec& x, double& err_cons, double& err_prec) {
    Vec grad;
    value = obj_fcn_alm(this, x, grad);
    //std::cout << " [d]:" << x.segment(idx_d, dim_d).transpose() << std::endl;
    double x_inf_norm = std::max(1.0, x.cwiseAbs().maxCoeff());
    double grad_inf_norm = grad.cwiseAbs().maxCoeff();
    double err_cone = 0;
    Vec3 prj;
    int n_cones = As.size();
    for (int k = 0; k < n_cones; ++k) {
      get_conic3_prj((1.0 / pa2.alm_rho) * mius[k] - As[k] * x - bs[k], prj);
      double err_cone1 =
          ((1.0 / pa2.alm_rho) * mius[k] - prj).cwiseAbs().maxCoeff();
      if (err_cone1 > err_cone) err_cone = err_cone1;
    }
    double err_ieq = ((-1.0 / pa2.alm_rho) * eta)
                         .cwiseMax(Aieq * x - bieq)
                         .cwiseAbs()
                         .maxCoeff();
    double err_eq = (G * x - h).cwiseAbs().maxCoeff();

    err_cons = std::max(std::max(err_cons, err_ieq), err_eq) / x_inf_norm;
    err_prec = grad_inf_norm / x_inf_norm;
    return value;
  }

  // ALM solving process
  int alm_solve(Vec& x1, double& value1, alm_param& pa,
                lbfgs::lbfgs_parameter_t pl) {
    pa2 = pa;
    pl2 = pl;
    int n_cones = As.size();
    for (int k = 0; k < n_cones; ++k) {
      mius[k].setZero();
    }
    eta.setZero();
    lambda.setZero();
    Vec x = x1;
    int lbfgs_ret, valid = 1;
    int alm_count = 0;
    double value;
    double err_cons, err_prec;
    value = get_alm_error(x, err_cons, err_prec);
    while (err_cons > pa2.alm_epsilon_cons || err_prec > pa2.alm_epsilon_prec) {
      // lbfgs inner layer optimize
      pl2.g_epsilon =
          std::max(std::pow(pa2.alm_xi, alm_count), pa2.alm_xi_min) *
          std::min(1.0, err_cons);
      pl2.g_epsilon = std::max(4e-5, pl2.g_epsilon);
      lbfgs_ret =
          lbfgs::lbfgs_optimize(x, value, obj_fcn_alm, nullptr, this, pl2);
      // alm value update
      int n_cones = As.size();
      for (int k = 0; k < n_cones; ++k) {
        get_conic3_prj(mius[k] - pa2.alm_rho * (As[k] * x + bs[k]), mius[k]);
      }
      eta = (eta + pa2.alm_rho * (Aieq * x - bieq)).cwiseMax(zeros_ieq);
      lambda += pa2.alm_rho * (G * x - h);
      // conditions udate
      value = get_alm_error(x, err_cons, err_prec);
      if (pa2.alm_iter_display > 0) {
        std::cout << "ALM iter: " << alm_count << " rho:" << pa2.alm_rho
                  << " lbfgs:" << lbfgs_ret << " ALM obj:" << value
                  << " err_cons:" << err_cons << " err_prec:" << err_prec
                  << std::endl;
      }
      if ((lbfgs_ret != 0 && lbfgs_ret != 1) || alm_count > pa2.alm_iter_max) {
        valid = 0;
        break;
      }
      // alm rho update
      ++alm_count;
      pa2.alm_rho = std::min(pa2.alm_beta, (1.0 + pa2.alm_gamma) * pa2.alm_rho);
    }
    x1 = x;
    value1 = value;
    return !valid;
  }

  // solve interface
  int solve(double& value, Vec& xa, Vec& xb, Vec& xc, Vec& xd) {
    // x0.segment(idx_a, dim_a) = Vec::Zero(dim_a);
    // x0.segment(idx_b, dim_b) = Vec::Ones(dim_b);
    // x0.segment(idx_c, dim_c) = Vec::Ones(dim_c);
    // x0.segment(idx_d, dim_d) = Vec::Ones(dim_d) * 0.5;
    x0.setZero();
    conicALMTOPP2::alm_param pa;
    lbfgs::lbfgs_parameter_t pl;
    int topp_ret = alm_solve(x0, value, pa, pl);
    get_result(xa, xb, xc, xd);
    value = c.dot(x0);
    return topp_ret;
  }

  // get result 
  int get_result(Vec& xa, Vec& xb, Vec& xc, Vec& xd){
    xa = x0.segment(idx_a, dim_a);
    xb = x0.segment(idx_b, dim_b);
    xc = x0.segment(idx_c, dim_c);
    xd = x0.segment(idx_d, dim_d);
    return 0;
  }
};

#endif  // CONIC_ALM_TOPP_01