#include <iomanip>
#include "SpaceHub/src/dev-tools.hpp"
#include "SpaceHub/src/macros.hpp"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/orbits/orbits.hpp"
#include "SpaceHub/src/stellar/stellar.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"

using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};

inline double critical_velocity_of_1_plus_2(double m1, double m2, double m3, double a) {
  double m_in = m1 + m2;
  double mu_out = m_in * m3 / (m_in + m3);

  return sqrt(space::consts::G / mu_out * m1 * m2 / a);
}

inline double critical_velocity_of_2_plus_2(double m1, double m2, double m3, double m4, double a_in, double a_out) {
  double m_in = m1 + m2;
  double m_out = m3 + m4;

  double mu_out = m_in * m_out / (m_in + m_out);

  return sqrt(space::consts::G / mu_out * (m1 * m2 / a_in + m3 * m4 / a_out));
}
/*
double get_max_b(double u, double v_inf, double rp) {
  double v2 = v_inf * v_inf;
  double tmp = u / v2 + rp;
  return sqrt(tmp * tmp - u * u / (v2 * v2));
}*/

double get_b_max(double v_c, double v_inf, double a_max) { return a_max * (4 * v_c / v_inf + 3); }

auto create_incident_orbit(double u, double v_inf, double b, double w, double i, double phi, double start_r) {
  using namespace space;

  double a = -u / (v_inf * v_inf);

  double e = sqrt(1 + b * b / (a * a));

  double l = a * (1 - e * e);

  double theta = acos((l - start_r) / (e * start_r));

  double y = start_r * sin(theta);

  double x = start_r * cos(theta);

  double v_slope = (x + a * e) / y * b * b / (a * a);

  double theta_v = atan(v_slope);

  double v = sqrt(u * (2 / start_r - 1 / a));

  Vector3d incid_pos{x, y * cos(phi), y * sin(phi)};

  Vector3d incid_vel{v * cos(theta_v), v * sin(theta_v) * cos(phi), v * sin(theta_v) * sin(phi)};

  orbit::euler_rotate(incid_pos, w, i, 0.0);

  orbit::euler_rotate(incid_vel, w, i, 0.0);

  return std::make_tuple(incid_pos, incid_vel);
}

double time_to_pericenter(double u, double v_inf, double b, double start_r) {
  double a = -u / (v_inf * v_inf);

  double e = sqrt(1 + b * b / (a * a));

  double l = a * (1 - e * e);

  double theta = acos((l - start_r) / (e * start_r));

  double cos_theta = cos(theta);

  double cosh_F = (e + cos_theta) / (1 + e * cos_theta);

  double F = acosh(cosh_F);

  double M = e * sinh(F) - F;

  double time_to_pericenter = sqrt(-a * a * a / u) * M;

  return time_to_pericenter;
}

template <typename Gen>
auto create_incident_M_star(Gen &gen, double u, double v_inf, double b, double start_r) {
  using namespace space;
  using namespace space::orbit;

  double w = randomGen::Uniform<double>::get(gen, 0, 2 * consts::pi);

  double i = acos(randomGen::Uniform<double>::get(gen, -1, 1));

  double phi = randomGen::Uniform<double>::get(gen, 0, 2 * consts::pi);

  auto [incid_pos, incid_vel] = create_incident_orbit(u, v_inf, b, w, i, phi, start_r);

  return std::make_tuple(incid_pos, incid_vel, w, i, phi);
}

double get_duration(double D, double v) { return D / v * 4; }

//------------------------------------fixed parameter-----------------------------------------------

template <typename Gen>
auto create_coplanner_incident_star(Gen &gen, double v, double b, double D) {
  using namespace space::orbit;
  using namespace space::unit;
  using namespace space::consts;

  Particle star{2 * m_solar, pow(2, 1.0 / 3) * r_solar};

  auto [incid_pos, incid_vel] = create_coplanner_orbit(gen, v, b, D, 0.0, 0.0, 0.0);

  move_particles_to(incid_pos, incid_vel, star);

  return std::make_tuple(star, v, b, 0.0, 0.0, 0.0);
}