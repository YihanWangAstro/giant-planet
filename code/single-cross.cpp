#include "scattering.hpp"

double V_DISPER{5 * space::unit::kms};
double AJ{5};
constexpr double DELTA{1e-5};

using namespace space::orbit;

void mono_single(std::string workdir, size_t idx, size_t sim_num, double m_dwarf) {
  using namespace space;
  std::random_device rd;
  std::mt19937 local_thread_gen(rd());

  std::fstream out_file{workdir + "_" + std::to_string(idx) + ".txt", std::fstream::out};

  // double v_inf = V_DISPER;  // space::randomGen::Maxwell<double>::get(local_thread_gen, V_DISPER);

  double a_j = AJ * space::unit::au;

  double v_inf = V_DISPER * space::unit::kms;

  double m_in = unit::m_solar + unit::m_jupiter;

  double mu_in = (unit::m_solar * unit::m_jupiter) / m_in;

  double u_out = space::consts::G * (m_dwarf + m_in);

  double v_c = critical_velocity_of_1_plus_2(unit::m_solar, unit::m_jupiter, m_dwarf, a_j);

  double b_max = get_b_max(v_c, v_inf, a_j);  /// get_max_b(u_out, V_DISPER / 5, 5 * a_j);

  double start_r = a_j * pow(2 * m_dwarf / (DELTA * mu_in), 1.0 / 3);

  for (size_t i = 0; i < sim_num; i++) {
    Particle sun{unit::m_solar, unit::r_solar}, jupiter{unit::m_jupiter, unit::r_jupiter};

    auto jupiter_orbit =
        Kepler{unit::m_solar, unit::m_jupiter, semi_latus_rectum(a_j, 0.0), 0.0, 0.0, 0.0, ISOTHERMAL, 0.0};

    move_particles_to(jupiter_orbit, jupiter);

    double b = sqrt(randomGen::Uniform<double>::get(local_thread_gen, 0, b_max * b_max));

    auto [incid_pos, incid_vel, w, incl, phi] = create_incident_M_star(local_thread_gen, u_out, v_inf, b, start_r);

    Particle star{m_dwarf, pow(m_dwarf, 1.0 / 3) * unit::r_solar};

    move_particles_to(incid_pos, incid_vel, star);

    move_to_com_coord(sun, jupiter);

    move_to_com_coord(sun, jupiter, star);

    double end_time = 0;

    if (v_inf > v_c) {
      end_time = 2 * time_to_pericenter(u_out, v_inf, b, start_r);
    } else {
      end_time = 20 * time_to_pericenter(u_out, v_inf, b, start_r);
    }

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    args.add_stop_point_operation(
        [&](auto &ptc) { space::display(out_file, i, jupiter_orbit, v_inf, b, w, incl, phi, ptc, "\r\n"); });

    spacex::SpaceXsim simulator{0, sun, jupiter, star};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  std::string output_name;

  size_t sim_num;

  space::tools::read_command_line(argc, argv, sim_num, V_DISPER, output_name);

  size_t thread_num = (std::thread::hardware_concurrency() > 1) ? std::thread::hardware_concurrency() : 1;

  /*double m_dwarf_min = 0.08 * space::unit::m_solar;

  double m_dwarf_max = 1 * space::unit::m_solar;

  double dm = (m_dwarf_max - m_dwarf_min) / thread_num;*/

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