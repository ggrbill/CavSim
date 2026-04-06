#include "Structures.hpp"


PrimeCoefficients::PrimeCoefficients(int n_x, int n_y)
    : alpha_x(n_x - 1, std::vector<Faces>(n_y)),
      beta_x(n_x - 1, std::vector<Faces>(n_y)),
      alpha_y(n_x, std::vector<Faces>(n_y - 1)),
      beta_y(n_x, std::vector<Faces>(n_y - 1)),
      Ap_u(n_x - 1, std::vector<double>(n_y, 0.0)),
      Aw_u(n_x - 1, std::vector<double>(n_y, 0.0)),
      Ae_u(n_x - 1, std::vector<double>(n_y, 0.0)),
      An_u(n_x - 1, std::vector<double>(n_y, 0.0)),
      As_u(n_x - 1, std::vector<double>(n_y, 0.0)),
      B_u(n_x - 1, std::vector<double>(n_y, 0.0)),
      Ap_v(n_x, std::vector<double>(n_y - 1, 0.0)),
      Aw_v(n_x, std::vector<double>(n_y - 1, 0.0)),
      Ae_v(n_x, std::vector<double>(n_y - 1, 0.0)),
      An_v(n_x, std::vector<double>(n_y - 1, 0.0)),
      As_v(n_x, std::vector<double>(n_y - 1, 0.0)),
      B_v(n_x, std::vector<double>(n_y - 1, 0.0)),
      Ap_p(n_x, std::vector<double>(n_y, 0.0)),
      Aw_p(n_x, std::vector<double>(n_y, 0.0)),
      Ae_p(n_x, std::vector<double>(n_y, 0.0)),
      An_p(n_x, std::vector<double>(n_y, 0.0)),
      As_p(n_x, std::vector<double>(n_y, 0.0)),
      B_p(n_x, std::vector<double>(n_y, 0.0)) {}

CavSimResult::CavSimResult(int n_x, int n_y)
    : u(n_x - 1, std::vector<double>(n_y, 0.0)),
      v(n_x, std::vector<double>(n_y - 1, 0.0)),
      u_old(n_x - 1, std::vector<double>(n_y, 0.0)),
      v_old(n_x, std::vector<double>(n_y - 1, 0.0)),
      u_hat(n_x - 1, std::vector<double>(n_y, 0.0)),
      v_hat(n_x, std::vector<double>(n_y - 1, 0.0)),
      P(n_x, std::vector<double>(n_y, 0.0)),
      Pn(n_x, std::vector<double>(n_y, 0.0)) {}
