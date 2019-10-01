#include "scattering.hpp"

double V_DISPER{30 * space::unit::kms};
double V_INF{5};
// double M_DWARF{0.2};
double AJ{1};

double AS{1};

constexpr double DELTA{1e-4};

using namespace space::orbit;

void mono_binary(std::string workdir, size_t idx, size_t sim_num, double m_dwarf) {
  using namespace space;
  std::random_device rd;
  std::mt19937 local_thread_gen(rd());

  std::fstream out_file{workdir + "_" + std::to_string(idx) + ".txt", std::fstream::out};
  // randomGen::PowerLaw<double>::get(local_thread_gen, -1.3, 0.08, 0.45);
  // double m_dwarf = M_DWARF * space::unit::m_solar;

  // double v_inf = V_INF * space::unit::kms;  // V_RATIO * critical_velocity_of_1_plus_2(sun.mass, jupiter.mass, );

  double a_j = AJ * space::unit::au;

  double a_s = AS * space::unit::au;

  double m_in = unit::m_solar + unit::m_jupiter;

  double mu_in = (unit::m_solar * unit::m_jupiter) / m_in;

  double u_out = space::consts::G * (m_dwarf + m_in);

  double b_max = 100 * space::unit::au;  // get_max_b(u_out, v_inf, 5 * a_j);

  double start_r = a_j * pow(2 * m_dwarf / (DELTA * mu_in), 1.0 / 3);

  for (size_t i = 0; i < sim_num; i++) {
    double v_inf = space::randomGen::Maxwell<double>::get(local_thread_gen, V_DISPER);

    Particle sun{unit::m_solar, unit::r_solar}, jupiter{unit::m_jupiter, unit::r_jupiter};

    auto jupiter_orbit =
        Kepler{unit::m_solar, unit::m_jupiter, semi_latus_rectum(a_j, 0.0), 0.0, 0.0, 0.0, ISOTHERMAL, 0.0};

    move_particles_to(jupiter_orbit, jupiter);

    move_to_com_coord(sun, jupiter);

    double b = sqrt(randomGen::Uniform<double>::get(local_thread_gen, 0, b_max * b_max));

    auto [incid_pos, incid_vel, w, incl, phi] = create_incident_M_star(local_thread_gen, u_out, v_inf, b, start_r);

    auto binary_orbit =
        Kepler{0.5 * m_dwarf, 0.5 * m_dwarf, semi_latus_rectum(a_s, 0.0), 0.0, 0.0, 0.0, ISOTHERMAL, 0.0};

    Particle star1{0.5 * m_dwarf, pow(0.5 * m_dwarf, 1.0 / 3) * unit::r_solar};

    Particle star2{0.5 * m_dwarf, pow(0.5 * m_dwarf, 1.0 / 3) * unit::r_solar};

    move_particles_to(binary_orbit, star2);

    move_particles_to(incid_pos, incid_vel, star1, star2);

    move_to_com_coord(sun, jupiter, star1, star2);

    double end_time = 2 * time_to_pericenter(u_out, v_inf, b, start_r);

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    args.add_stop_point_operation([&](auto &ptc) {
      // out_file << PACK(i, ' ', jupiter_orbit, ' ', v_inf, ' ', b, ' ', w, ' ', incl, ' ', phi, ' ', ptc, ' ',
      // "\r\n");
      space::display(out_file, i, jupiter_orbit, v_inf, b, w, incl, phi, ptc, binary_orbit, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  std::string output_name;

  size_t sim_num;

  space::tools::read_command_line(argc, argv, sim_num, AS, AJ, output_name);

  size_t thread_num = 80;

  double m_dwarf_min = 0.08 * space::unit::m_solar;

  double m_dwarf_max = 1 * space::unit::m_solar;

  double dm = (m_dwarf_max - m_dwarf_min) / thread_num;

  double m_dwarf = m_dwarf_min;

  std::vector<std::thread> threads;

  for (size_t i = 0; i < thread_num; ++i) {
    threads.emplace_back(std::thread{mono_single, output_name, i, sim_num, m_dwarf});
    m_dwarf += dm;
  }

  for (auto &th : threads) {
    th.join();
  }
  // space::multiThread::auto_multi_thread(mono_single, out_file, sim_num);

  return 0;
}