#include <iomanip>
#include "SpaceHub/src/dev-tools.hpp"
#include "SpaceHub/src/macros.hpp"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/orbits/orbits.hpp"
#include "SpaceHub/src/stellar/stellar.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"

using ConFile = space::multiThread::ConcurrentFile;
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

auto create_1_planet_system(double aj, double ej = 0) {
  using namespace space::orbit;
  using namespace space;

  Particle sun{unit::m_solar, unit::r_solar}, jupiter{unit::m_jupiter, unit::r_jupiter};

  auto jupiter_orbit = Kepler{unit::m_solar, unit::m_jupiter, semi_latus_rectum(aj, ej), ej, 0.0, 0.0, ISOTHERMAL, 0.0};

  move_particles_to(jupiter_orbit, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}

template <typename Gen>
auto create_incident_orbit(Gen &gen, double v_inf, double b, double D, double w, double i, double phi) {
  using namespace space;

  Vector3d incid_pos{D, b * cos(phi), b * sin(phi)};

  Vector3d incid_vel{-v_inf, 0, 0.0};

  orbit::euler_rotate(incid_pos, w, i, 0.0);

  orbit::euler_rotate(incid_vel, w, i, 0.0);

  return std::make_tuple(incid_pos, incid_vel);
}

template <typename Gen>
auto create_incident_M_star(Gen &gen, double v_dispersion, double b_max, double D) {
  using namespace space;
  using namespace space::orbit;

  double m_star = randomGen::PowerLaw<double>::get(gen, -1.3, 0.08, 0.45);

  Particle star{m_star, pow(m_star, 1.0 / 3) * unit::r_solar};

  double v_inf = randomGen::Maxwell<double>::get(gen, v_dispersion);

  double b = sqrt(randomGen::Uniform<double>::get(gen, 0, b_max * b_max));

  double w = randomGen::Uniform<double>::get(gen, 0, 2 * consts::pi);

  double i = acos(randomGen::Uniform<double>::get(gen, -1, 1));

  double phi = randomGen::Uniform<double>::get(gen, 0, 2 * consts::pi);

  auto [incid_pos, incid_vel] = create_incident_orbit(gen, v_inf, b, D, w, i, phi);

  move_particles_to(incid_pos, incid_vel, star);

  return std::make_tuple(star, incid_vel.norm(), b, w, i, phi);
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