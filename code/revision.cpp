#include "scattering.hpp"

double V_MEAN{5};
double V_DISPER{5};
double B_MAX{100};
double B_MIN{-100};
double AJ{1};
using namespace space::orbit;

void mono_single(ConFile out_file, size_t sim_num) {
  std::random_device rd;
  std::mt19937 local_thread_gen(rd());
  double D = 500 * space::unit::au;

  for (size_t i = 0; i < sim_num; i++) {
    double b_var = space::randomGen::Uniform<double>::get(B_MIN, B_MAX);

    auto [sun, jupiter, jupiter_orbit] = create_1_planet_system(AJ, 0);

    auto [star, v_inf, b, w, incl, phi] = create_incident_M_star(local_thread_gen, V_MEAN, b_var, D);

    move_to_com_coord(sun, jupiter);
    move_to_com_coord(sun, jupiter, star);

    double end_time = get_duration(D, v_inf);

    spacex::SpaceXsim::RunArgs args;

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    args.add_stop_point_operation([&](auto &ptc) {
      out_file << PACK(i, ' ', jupiter_orbit, ' ', v_inf, ' ', b, ' ', w, ' ', incl, ' ', phi, ' ', ptc, ' ', "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star};

    simulator.run(args);
  }
}

int main(int argc, char **argv) {
  std::string output_name;

  size_t sim_num;

  space::tools::read_command_line(argc, argv, sim_num, B_MIN, B_MAX, AJ, output_name);

  V_MEAN *= space::unit::kms;
  V_DISPER = V_MEAN * sqrt(space::consts::pi / 8.0);
  B_MIN *= space::unit::au;
  B_MAX *= space::unit::au;
  AJ *= space::unit::au;

  std::string file_name = output_name + ".txt";

  auto out_file = space::multiThread::make_thread_safe_fstream(file_name, std::fstream::out);

  out_file << std::setprecision(6);

  space::multiThread::auto_multi_thread(mono_single, out_file, sim_num);

  return 0;
}