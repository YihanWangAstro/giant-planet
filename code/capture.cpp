#include "SpaceHub/src/macros.hpp"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/scattering/cross-section.hpp"
#include "SpaceHub/src/stellar/stellar.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceX/SpaceHubWrapper.hpp"

using ConFile = space::multi_thread::ConcurrentFile;
using Particle = spacex::SpaceXsim::Particle;
using ParticleSys = spacex::SpaceXsim::RunArgs::ParticleSys;

using namespace space;
using namespace space::orbit;

auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};

void mono_single(size_t thread_idx, std::string workdir, size_t sim_num, double a_j, double v_inf) {
  std::fstream out_file{workdir + "_" + std::to_string(thread_idx) + ".txt", std::fstream::out};
  double e_j = 0;
  double m_d = 0.2 * unit::m_solar;
  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_d);
  double const delta = 1e-5;

  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{unit::m_solar, unit::r_solar}, jupiter{unit::m_jupiter, unit::r_jupiter}, dwarf{m_d, r_d};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, e_j, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter), dwarf, v_inf, delta);

    move_particles(in_orbit, dwarf);

    move_to_COM_frame(sun, jupiter, dwarf);

    double end_time = ((E_tot(sun, jupiter, dwarf) > 0) ? 2.0 : 20.0) * time_to_periapsis(cluster(sun, jupiter), dwarf);

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { space::display(out_file, i, jupiter_orbit, v_inf, in_orbit, ptc, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, jupiter, dwarf};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  std::string output_name;

  size_t sim_num;

  double a_j{5};

  double v_inf{1};

  tools::read_command_line(argc, argv, sim_num, output_name, a_j, v_inf);

  a_j *= unit::au;

  v_inf *= unit::kms;

  multi_thread::auto_indexed_multi_thread(mono_single, output_name, sim_num, a_j, v_inf);

  return 0;
}