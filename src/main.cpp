#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Vec3.hpp"
#include "Particle.hpp"
#include "Fields.hpp"
#include "Integrator.hpp"

int main()
{
    std::cout << "Charged Particle Simulator\n";

    // --------------------------------------------------
    // Normalized simulation parameters
    //
    // Reference normalization:
    // Omega_0 = |q| B_0 / m
    // t_0     = 1 / Omega_0
    // rho_0   = v_0 / Omega_0
    //
    // Therefore, for this benchmark:
    // m = 1, |q| = 1, B_0 = 1
    // --------------------------------------------------

    Particle particle{
        {0.0, 0.0, 0.0},   // position
        {1.0, 0.0, 0.2},   // velocity
        1.0,                // normalized mass
        1.0                 // normalized charge
    };

    FieldConfig fields;
    fields.type = FieldCase::UniformB;
    fields.B0 = 1.0;

    const double dt = 0.01;
    const double final_time = 20.0;

    if (particle.mass <= 0.0) {
        std::cerr << "Error: particle mass must be positive.\n";
        return 1;
    }

    if (dt <= 0.0) {
        std::cerr << "Error: time step must be positive.\n";
        return 1;
    }

    if (fields.B0 == 0.0) {
        std::cerr << "Error: uniform-B benchmark requires B0 != 0.\n";
        return 1;
    }

    std::ofstream output("results/uniform_B.csv");

    if (!output) {
        std::cerr << "Error: could not open output file.\n";
        return 1;
    }

    if (particle.charge == 0.0) {
        std::cerr << "Error: uniform-B benchmark requires nonzero charge.\n";
        return 1;
    }

    output << std::setprecision(12);

    output
        << "t,"
        << "x,y,z,"
        << "vx,vy,vz,"
        << "kinetic_energy\n";

    double time = 0.0;

    const double pi = std::acos(-1.0);

    const double initial_energy =
        particle.kinetic_energy();

    const double gyrofrequency =
        std::abs(particle.charge)
        * std::abs(fields.B0)
        / particle.mass;

    const double gyroperiod =
        2.0 * pi / gyrofrequency;

    const double perpendicular_speed =
        std::sqrt(
            particle.velocity.x * particle.velocity.x
            + particle.velocity.y * particle.velocity.y
        );

    const double larmor_radius =
        perpendicular_speed / gyrofrequency;

    double max_relative_energy_error = 0.0;

    while (time <= final_time) {

        const double current_energy =
            particle.kinetic_energy();

        const double relative_energy_error =
            std::abs(current_energy - initial_energy)
            / initial_energy;

        if (relative_energy_error > max_relative_energy_error) {
            max_relative_energy_error = relative_energy_error;
        }

        output
            << time << ","
            << particle.position.x << ","
            << particle.position.y << ","
            << particle.position.z << ","
            << particle.velocity.x << ","
            << particle.velocity.y << ","
            << particle.velocity.z << ","
            << particle.kinetic_energy()
            << '\n';

        rk4_step(
            particle,
            fields,
            time,
            dt
        );

        time += dt;
    }
    
    output.close();

    const double final_energy =
    particle.kinetic_energy();

    std::cout
        << "\n--- Uniform-B benchmark ---\n"
        << "Gyrofrequency:             "
        << gyrofrequency << '\n'
        << "Gyroperiod:                "
        << gyroperiod << '\n'
        << "Larmor radius:             "
        << larmor_radius << '\n'
        << "Initial kinetic energy:    "
        << initial_energy << '\n'
        << "Final kinetic energy:      "
        << final_energy << '\n'
        << "Maximum relative E error:  "
        << max_relative_energy_error << '\n';

    std::cout
        << "Simulation finished.\n"
        << "Output: results/uniform_B.csv\n";

    return 0;
}
