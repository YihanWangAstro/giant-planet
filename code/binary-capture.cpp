#include "SpaceHub/src/args-callback/callbacks.hpp"
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
using namespace unit;
auto collision = [](auto &ptc) -> bool {
  size_t number = ptc.number();
  for (size_t i = 0; i < number; ++i) {
    for (size_t j = i + 1; j < number; ++j) {
      if (ptc.Collide(i, j)) return true;
    }
  }
  return false;
};

void mono_binary(size_t thread_idx, std::string workdir, size_t sim_num, double a_j, double v_inf, double a_s) {
  std::fstream out_file{workdir + "_" + std::to_string(thread_idx) + ".txt", std::fstream::out};
  double e_j = 0;
  double m_d = 0.2_Ms;
  double r_d = stellar::stellar_radius(stellar::StarType::STAR, m_d);
  double const delta = 1e-5;
  double e_s = 0;
  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun1{1_Ms, 1_Rs}, sun2{1_Ms, 1_Rs}, jupiter{1_Mj, 1_Rj}, dwarf{m_d, r_d};

    auto jupiter_orbit = EllipOrbit{sun1.mass, jupiter.mass, a_j, e_j, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    auto solar_orbit = EllipOrbit{sun1.mass, sun2.mass, a_s, e_s, isotherm, isotherm, isotherm, isotherm};

    move_particles(solar_orbit, sun1, jupiter);

    move_to_COM_frame(sun1, sun2, jupiter);

    auto in_orbit = scattering::incident_orbit(cluster(sun1, sun2, jupiter), dwarf, v_inf, delta);

    move_particles(in_orbit, dwarf);

    move_to_COM_frame(sun1, sun2, jupiter, dwarf);

    double end_time =
        ((E_tot(sun1, sun2, jupiter, dwarf) > 0) ? 2.0 : 20.0) * time_to_periapsis(cluster(sun1, sun2, jupiter), dwarf);

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { space::display(out_file, i, jupiter_orbit, v_inf, in_orbit, ptc, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun1, sun2, jupiter, dwarf};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  size_t sim_num;
  std::string output_name;
  double a_j, v_inf, a_s;

  tools::read_command_line(argc, argv, sim_num, output_name, a_j, v_inf, a_s);

  a_j *= unit::AU;
  a_s *= unit::AU;
  v_inf *= unit::kms;

  auto th_num = multi_thread::machine_thread_num;

  multi_thread::indexed_multi_thread(th_num, mono_binary, output_name, sim_num, a_j, v_inf, a_s);

  return 0;
}