#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "Fields.hpp"
#include "Integrator.hpp"
#include "Particle.hpp"

int main()
{
    std::cout << "Charged Particle Simulator\n";

    // --------------------------------------------------
    // Benchmark configuration
    // --------------------------------------------------
    // All quantities are expressed in normalized units.
    // The active setup below corresponds to the nonuniform-E benchmark.
    // Parameter sets for the other benchmarks are documented in the report.

    Particle particle{
        {0.0, 0.0, 0.0},   // Initial position
        {1.0, 0.0, 0.2},   // Initial velocity
        1.0,                // Normalized mass
        1.0                 // Normalized charge
    };

    FieldConfig fields;
    fields.type = FieldCase::NonuniformE;
    fields.B0 = 1.0;
    fields.E0 = 0.02;
    fields.wavenumber = 0.4;

    const double dt = 0.01;
    const double final_time = 60.0;

    const std::string output_path =
        "results/nonuniform_E.csv";

    // --------------------------------------------------
    // Basic input validation
    // --------------------------------------------------

    if (particle.mass <= 0.0) {
        std::cerr << "Error: particle mass must be positive.\n";
        return 1;
    }

    if (particle.charge == 0.0) {
        std::cerr << "Error: particle charge must be nonzero.\n";
        return 1;
    }

    if (dt <= 0.0) {
        std::cerr << "Error: time step must be positive.\n";
        return 1;
    }

    if (final_time < 0.0) {
        std::cerr << "Error: final time must be nonnegative.\n";
        return 1;
    }

    // --------------------------------------------------
    // Output file
    // --------------------------------------------------

    std::filesystem::create_directories("results");

    std::ofstream output(output_path);

    if (!output) {
        std::cerr << "Error: could not open output file.\n";
        return 1;
    }

    output << std::setprecision(12);
    output << "t,x,y,z,vx,vy,vz,kinetic_energy\n";

    // --------------------------------------------------
    // Time integration
    // --------------------------------------------------

    const std::size_t num_steps =
        static_cast<std::size_t>(
            std::llround(final_time / dt)
        );

    for (std::size_t step = 0; step <= num_steps; ++step) {
        const double time =
            static_cast<double>(step) * dt;

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

        if (step < num_steps) {
            rk4_step(
                particle,
                fields,
                time,
                dt
            );
        }
    }

    std::cout
        << "Simulation finished.\n"
        << "Output: " << output_path << '\n';

    return 0;
}
