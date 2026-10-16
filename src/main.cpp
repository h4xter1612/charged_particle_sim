#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "Vec3.hpp"
#include "Particle.hpp"
#include "Fields.hpp"
#include "Integrator.hpp"

int main()
{
    std::cout << "Charged Particle Simulator\n";

    // --------------------------------------------------
    // Particle initial conditions
    // --------------------------------------------------

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

    // --------------------------------------------------
    // Field configuration
    // --------------------------------------------------

    FieldConfig fields;

    fields.type = FieldCase::ExB;
    fields.B0 = 1.0;
    fields.E0 = 0.1;

    // --------------------------------------------------
    // Simulation parameters
    // --------------------------------------------------

    const double dt = 0.01;
    const double final_time = 20.0;

    const std::string output_path =
        "results/exb.csv";

    // --------------------------------------------------
    // Basic input validation
    // --------------------------------------------------

    if (particle.mass <= 0.0) {
        std::cerr
            << "Error: particle mass must be positive.\n";
        return 1;
    }

    if (particle.charge == 0.0) {
        std::cerr
            << "Error: particle charge must be nonzero.\n";
        return 1;
    }

    if (dt <= 0.0) {
        std::cerr
            << "Error: time step must be positive.\n";
        return 1;
    }

    // --------------------------------------------------
    // Output file
    // --------------------------------------------------

    std::ofstream output(output_path);

    if (!output) {
        std::cerr
            << "Error: could not open output file.\n";
        return 1;
    }

    output << std::setprecision(12);

    output
        << "t,"
        << "x,y,z,"
        << "vx,vy,vz,"
        << "kinetic_energy\n";

    // --------------------------------------------------
    // Time integration
    // --------------------------------------------------

    double time = 0.0;

    while (time <= final_time) {

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

    std::cout
        << "Simulation finished.\n"
        << "Output: " << output_path << '\n';

    return 0;
}
