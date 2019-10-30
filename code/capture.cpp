#include "SpaceHub/src/orbits/orbits.hpp"
#include "SpaceHub/src/orbits/particle-manip.hpp"
#include "SpaceHub/src/scattering/cross-section.hpp"

#include "SpaceHub/src/dev-tools.hpp"
#include "SpaceHub/src/macros.hpp"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"

#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
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

void mono_single(std::string workdir, size_t idx, size_t sim_num, double m_dwarf) {
  std::fstream out_file{workdir + "_" + std::to_string(idx) + ".txt", std::fstream::out};

  double a_j = 1;
  double e_j = 0;
  double m_d = 0.2;
  double r_d = 1;
  double v_inf = 1;
  double delta = 1e-5;
  for (size_t i = 0; i < sim_num; ++i) {
    Particle sun{unit::m_solar, unit::r_solar}, jupiter{unit::m_jupiter, unit::r_jupiter}, m_dwarf{m_d, r_d};

    auto jupiter_orbit = EllipOrbit{sun.mass, jupiter.mass, a_j, e_j, isotherm, isotherm, isotherm, isotherm};

    move_particles(jupiter_orbit, jupiter);

    auto in_orbit = scattering::incident_orbit(cluster(sun, jupiter), m_dwarf, v_inf, 1e-5);

    move_particles(in_orbit, m_dwarf);

    move_to_COM_frame(sun, jupiter, m_dwarf);

    double end_time = time_to_periapsis(cluster(sun, jupiter), m_dwarf);

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { space::display(out_file, i, jupiter_orbit, v_inf, in_orbit, ptc, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, jupiter, m_dwarf};

    simulator.run(args);
  }

  int main(int argc, char **argv) {
    Particle sun{unit::m_solar, unit::r_solar}, jupiter{unit::m_jupiter, unit::r_jupiter};

    auto x = cluster(sun, jupiter);

    if constexpr (is_ranges_v<decltype(x)>) {
      std::cout << "is container";
    }
    std::cout << x;

    std::string output_name;

    size_t sim_num;

    space::tools::read_command_line(argc, argv, sim_num, V_DISPER, output_name);

    size_t thread_num = (std::thread::hardware_concurrency() > 1) ? std::thread::hardware_concurrency() : 1;

    double m_dwarf = 0.2 * space::unit::m_solar;

    std::vector<std::thread> threads;

    for (size_t i = 0; i < thread_num; ++i) {
      threads.emplace_back(std::thread{mono_single, output_name, i, sim_num, m_dwarf});
    }

    for (auto &th : threads) {
      th.join();
    }
    return 0;
  }