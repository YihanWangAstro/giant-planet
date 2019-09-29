#include "SpaceHub/src/dev-tools.hpp"
#include "SpaceHub/src/macros.hpp"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/orbits/orbits.hpp"
#include "SpaceHub/src/stellar/stellar.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"
#include <iomanip>

using ConFile = space::multiThread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

using namespace space::orbit;

struct Coord {
  Vector3d pos, vel;
  friend std::ostream &operator<<(std::ostream &os, Coord const &c) {
    os << c.pos.x << ' ' << c.pos.y << ' ' << c.pos.z << ' ' << c.vel.x << ' '
       << c.vel.y << ' ' << c.vel.z;
    return os;
  }
};

template <size_t num> struct Coords_ {
  std::array<Coord, num> data;

  // template <size_t S>
  friend std::ostream &operator<<(std::ostream &os, Coords_<num> const &d) {
    for (auto const &c : d.data) {
      os << c << ' ';
    }
    return os;
  }
};

auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j))
        return true;
    }
  }
  return false;
};

template <typename Coords, typename PTCS>
void save_states(PTCS const &ptc, Coords &c) {
  size_t num = c.data.size();

  for (size_t i = 0; i < num; ++i) {
    c.data[i].pos = ptc[i].pos;
    c.data[i].vel = ptc[i].vel;
  }
}

auto create_1_planet_system(double aj, double ej = 0) {
  using namespace space::orbit;
  using namespace space::unit;
  using namespace space::consts;
  Particle sun{m_solar, r_solar}, jupiter{m_jupiter, r_jupiter};

  auto jupiter_orbit =
      Kepler{m_solar,    m_jupiter, semi_latus_rectum(aj, ej), ej, 0.0, 0.0,
             ISOTHERMAL, 0.0};

  move_particles_to(jupiter_orbit, jupiter);

  return std::make_tuple(sun, jupiter, jupiter_orbit);
}

auto create_2_planet_system(double aj, double an, double ej = 0,
                            double en = 0) {
  using namespace space::orbit;
  using namespace space::unit;
  using namespace space::consts;
  Particle sun{m_solar, r_solar}, jupiter{m_jupiter, r_jupiter},
      neptune{m_neptune, r_neptune};

  Kepler jupiter_orbit, neptune_orbit;

  jupiter_orbit =
      Kepler{m_solar,    m_jupiter, semi_latus_rectum(aj, ej), ej, 0.0, 0.0,
             ISOTHERMAL, 0.0};

  neptune_orbit =
      Kepler{m_solar,    m_neptune, semi_latus_rectum(an, en), en, 0.0, 0.0,
             ISOTHERMAL, 0.0};

  move_particles_to(jupiter_orbit, jupiter);

  move_particles_to(neptune_orbit, neptune);

  return std::make_tuple(sun, jupiter, neptune, jupiter_orbit, neptune_orbit);
}

/*
template<typename Gen>
auto create_incident_orbit(Gen &gen, double v_dispersion, double b_max, double
D){ double v_inf = space::randomGen::Maxwell<double>::get(gen, v_dispersion);

    double b = space::randomGen::Uniform<double>::get(gen, 0, b_max);

    double w = space::randomGen::Uniform<double>::get(gen, 0,
2*space::consts::pi);

    double i = acos(space::randomGen::Uniform<double>::get(gen, -1, 1));;

    double phi = space::randomGen::Uniform<double>::get(gen, 0,
2*space::consts::pi);

    Vector3d incid_pos{D, b*cos(phi), b*sin(phi)};

    Vector3d incid_vel{-v_inf, 0, 0.0};

    euler_rotate(incid_pos, w, i, 0.0);

    euler_rotate(incid_vel, w, i, 0.0);

    return std::make_tuple(incid_pos, incid_vel, b, w, i, phi);
}

template<typename Gen>
auto create_incident_star(Gen &gen, double v_dispersion, double b_max, double D)
{ using namespace space::orbit; using namespace space::unit; using namespace
space::consts;

    Particle star{2*m_solar, pow(2,1.0/3)*r_solar};

    auto [incid_pos, incid_vel, b, w, i, phi] = create_incident_orbit(gen,
v_dispersion, b_max, D);

    move_particles_to(incid_pos, incid_vel, star);

    return std::make_tuple(star, incid_vel.norm(), b, w, i, phi);
}

template<typename Gen>
auto create_incident_binary(Gen &gen, double v_dispersion, double b_max, double
D, double as = space::unit::au, double es = 0) { using namespace space::orbit;
    using namespace space::unit;
    using namespace space::consts;

    Particle star1{m_solar, r_solar}, star2{m_solar, r_solar};

    Kepler star_orbit = Kepler{m_solar, m_solar, semi_latus_rectum(as, es), es,
ISOTHERMAL, ISOTHERMAL, ISOTHERMAL, ISOTHERMAL};

    move_particles_to(star_orbit, star2);

    auto [incid_pos, incid_vel, b, w, i, phi] = create_incident_orbit(gen,
v_dispersion, b_max, D);

    move_particles_to(incid_pos, incid_vel, star1, star2);

    return std::make_tuple(star1, star2, star_orbit, incid_vel.norm(), b, w, i,
phi);
}*/

double get_duration(double D, double v) { return D / v * 4; }

//------------------------------------fixed
// parameter-----------------------------------------------
template <typename Gen>
auto create_coplanner_orbit(Gen &gen, double v, double b, double D) {
  return std::make_tuple(Vector3d{D, b, 0.0}, Vector3d{-v, 0.0, 0.0});
}

template <typename Gen>
auto create_coplanner_incident_binary(Gen &gen, double v, double b, double D,
                                      double as_min, double as_max, double I,
                                      double es) {
  using namespace space::orbit;
  using namespace space::unit;
  using namespace space::consts;

  Particle star1{m_solar, r_solar}, star2{m_solar, r_solar};

  double as = space::randomGen::Logarithm<double>::get(gen, as_min, as_max);

  Kepler star_orbit =
      Kepler{m_solar, m_solar, semi_latus_rectum(as, es), es, I, 0.0, 0.0, 0.0};

  move_particles_to(star_orbit, star2);

  auto [incid_pos, incid_vel] = create_coplanner_orbit(gen, v, b, D);

  move_particles_to(incid_pos, incid_vel, star1, star2);

  return std::make_tuple(star1, star2, star_orbit, v, b, 0.0, 0.0, 0.0);
}

template <typename Gen>
auto create_coplanner_incident_binary(Gen &gen, double v, double b, double D,
                                      double as, double I, double es) {
  using namespace space::orbit;
  using namespace space::unit;
  using namespace space::consts;

  Particle star1{m_solar, r_solar}, star2{m_solar, r_solar};

  Kepler star_orbit = Kepler{
      m_solar, m_solar, semi_latus_rectum(as, es), es, I, 0.0, 0.0, ISOTHERMAL};

  move_particles_to(star_orbit, star2);

  auto [incid_pos, incid_vel] = create_coplanner_orbit(gen, v, b, D);

  move_particles_to(incid_pos, incid_vel, star1, star2);

  return std::make_tuple(star1, star2, star_orbit, v, b, 0.0, 0.0, 0.0);
}

template <typename Gen>
auto create_coplanner_incident_star(Gen &gen, double v, double b, double D) {
  using namespace space::orbit;
  using namespace space::unit;
  using namespace space::consts;

  Particle star{2 * m_solar, pow(2, 1.0 / 3) * r_solar};

  auto [incid_pos, incid_vel] = create_coplanner_orbit(gen, v, b, D);

  move_particles_to(incid_pos, incid_vel, star);

  return std::make_tuple(star, v, b, 0.0, 0.0, 0.0);
}