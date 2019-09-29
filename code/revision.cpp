#include "scattering.hpp"

double V_MEAN{5};
double B_MAX{100};
double B_MIN{-100};
double AJ{1};
double AN{5};
double AS{1};

void mono_single(ConFile out_file, size_t sim_num) {
  std::random_device rd;
  std::mt19937 local_thread_gen(rd());
  double D = 500 * space::unit::au;

  for (size_t i = 0; i < sim_num; i++) {
    double b_var = space::randomGen::Uniform<double>::get(B_MIN, B_MAX);

    auto [sun, jupiter, jupiter_orbit] = create_1_planet_system(AJ, 0);

    Coords_<3> state_close;

    auto [star, v_inf, b, w, incl, phi] =
        create_coplanner_incident_star(local_thread_gen, V_MEAN, b_var, D);

    move_to_com_coord(sun, jupiter);
    move_to_com_coord(sun, jupiter, star);

    double end_time = get_duration(D, v_inf);

    spacex::SpaceXsim::RunArgs args;

    double closest_r{1000 * space::unit::pc};

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    double solar_tot_m = total_mass(sun, jupiter);

    args.add_post_step_operation([&](auto const &ptc) {
      auto CoM = (ptc[0].m * ptc[0].pos + ptc[1].m * ptc[1].pos) / solar_tot_m;

      auto dss = Distance(CoM, ptc[2].pos);
      if (dss < closest_r) {
        save_states(ptc, state_close);
      }
    });

    args.add_stop_point_operation([&](auto &ptc) {
      out_file << PACK(i, ' ', jupiter_orbit, ' ', v_inf, ' ', b, ' ', w, ' ',
                       incl, ' ', phi, ' ', ptc, ' ', state_close, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star};

    simulator.run(args);
  }
}

void multi_single(ConFile out_file, size_t sim_num) {
  std::random_device rd;
  std::mt19937 local_thread_gen(rd());
  double D = 500 * space::unit::au;

  for (size_t i = 0; i < sim_num; i++) {
    double b_var = space::randomGen::Uniform<double>::get(B_MIN, B_MAX);

    auto [sun, jupiter, neptune, jupiter_orbit, neptune_orbit] =
        create_2_planet_system(AJ, AN, 0, 0);

    Coords_<4> state_close;

    auto [star, v_inf, b, w, incl, phi] =
        create_coplanner_incident_star(local_thread_gen, V_MEAN, b_var, D);

    move_to_com_coord(sun, jupiter, neptune);
    move_to_com_coord(sun, jupiter, neptune, star);

    double end_time = get_duration(D, v_inf);

    spacex::SpaceXsim::RunArgs args;

    double closest_r{1000 * space::unit::pc};

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    double solar_tot_m = total_mass(sun, jupiter, neptune);

    args.add_post_step_operation([&](auto const &ptc) {
      auto CoM = (ptc[0].m * ptc[0].pos + ptc[1].m * ptc[1].pos +
                  ptc[2].m * ptc[2].pos) /
                 solar_tot_m;

      auto dss = Distance(CoM, ptc[3].pos);
      if (dss < closest_r) {
        save_states(ptc, state_close);
      }
    });

    args.add_stop_point_operation([&](auto &ptc) {
      out_file << PACK(i, ' ', jupiter_orbit, ' ', neptune_orbit, ' ', v_inf,
                       ' ', b, ' ', w, ' ', incl, ' ', phi, ' ', ptc, ' ',
                       state_close, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, neptune, star};

    simulator.run(args);
  }
}

void mono_binary(ConFile out_file, size_t sim_num) {
  std::random_device rd;
  std::mt19937 local_thread_gen(rd());
  double D = 500 * space::unit::au;

  double I = 0;

  for (size_t i = 0; i < sim_num; i++) {
    double b_var = space::randomGen::Uniform<double>::get(B_MIN, B_MAX);

    auto [sun, jupiter, jupiter_orbit] = create_1_planet_system(AJ, 0);

    Coords_<4> state_close;

    auto [star1, star2, star_orbit, v_inf, b, w, incl, phi] =
        create_coplanner_incident_binary(local_thread_gen, V_MEAN, b_var, D, AS,
                                         I, 0.0);

    move_to_com_coord(sun, jupiter);
    move_to_com_coord(sun, jupiter, star1, star2);

    double end_time = get_duration(D, v_inf);

    spacex::SpaceXsim::RunArgs args;

    double closest_r{1000 * space::unit::pc};

    args.add_stop_condition(collision);

    args.add_stop_condition(end_time);

    double solar_tot_m = total_mass(sun, jupiter);

    double bi_tot_m = total_mass(star1, star2);

    args.add_post_step_operation([&](auto const &ptc) {
      auto CoM = (ptc[0].m * ptc[0].pos + ptc[1].m * ptc[1].pos) / solar_tot_m;

      auto bCoM = (ptc[2].m * ptc[2].pos + ptc[3].m * ptc[3].pos) / bi_tot_m;

      auto dss = Distance(CoM, bCoM);
      if (dss < closest_r) {
        save_states(ptc, state_close);
      }
    });

    args.add_stop_point_operation([&](auto &ptc) {
      out_file << PACK(i, ' ', jupiter_orbit, ' ', star_orbit, ' ', v_inf, ' ',
                       b, ' ', w, ' ', incl, ' ', phi, ' ', ptc, ' ',
                       state_close, "\r\n");
    });

    spacex::SpaceXsim simulator{0, sun, jupiter, star1, star2};

    simulator.run(args);

    double random =
        space::randomGen::Uniform<double>::get(local_thread_gen, 0, 1);

    if (random <= 0.5) {
      I = 0;
    } else {
      I = space::consts::pi;
    }
  }
}

int main(int argc, char **argv) {
  std::string output_name;

  std::string sim_type;

  size_t sim_num;

  space::tools::read_command_line(argc, argv, sim_type, sim_num, B_MIN, B_MAX,
                                  AJ, AN, AS, output_name);

  V_MEAN *= space::unit::kms;
  B_MIN *= space::unit::au;
  B_MAX *= space::unit::au;
  AJ *= space::unit::au;
  AN *= space::unit::au;
  AS *= space::unit::au;

  std::string file_name = output_name + ".txt";

  auto out_file = space::multiThread::make_thread_safe_fstream(
      file_name, std::fstream::out);

  out_file << std::setprecision(6);

  if (sim_type == "mono-single") {
    space::multiThread::auto_multi_thread(mono_single, out_file, sim_num);
  } else if (sim_type == "multi-single") {
    space::multiThread::auto_multi_thread(multi_single, out_file, sim_num);
  } else if (sim_type == "mono-binary") {
    space::multiThread::auto_multi_thread(mono_binary, out_file, sim_num);
  } else {
    std::cout << "incorrect simulation type!\n";
  }

  return 0;
}