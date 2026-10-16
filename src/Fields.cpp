#include "Fields.hpp"
#include <stdexcept>
#include <cmath>

Vec3 electric_field(
    const Vec3& position,
    double time,
    const FieldConfig& config
)
{
    // Evita warnings por ahora.
    (void)position;
    (void)time;

    switch (config.type) {

        case FieldCase::UniformB:
            // E = 0
            return {0.0, 0.0, 0.0};


        case FieldCase::ExB:
            // Uniform electric field:
            // E = E0 x-hat
            return {config.E0, 0.0, 0.0};


        case FieldCase::GradB:
            // No electric field for the grad-B benchmark.
            return {0.0, 0.0, 0.0};


        case FieldCase::Curvature:
            // No electric field for the curvature-drift benchmark.
            return {0.0, 0.0, 0.0};


        case FieldCase::NonuniformE: {
            // Spatially nonuniform electric field:
            //
            // E(y) = E0 * cos(k x) x-hat
            //
            // For a finite Larmor radius, the particle samples different values of E during each gyro-orbit.

            const double Ex = config.E0 * std::cos(config.wavenumber * position.x);

            return {
                Ex,
                0.0,
                0.0
            };
        }

        case FieldCase::Polarization: {
            // Slowly time-varying perpendicular electric field:
            //
            // E(t) = E0 * cos(omega * t) x-hat
            //
            // This produces both a time-dependent E x B curvature-drift and a polarization drift.

            const double Ex = config.E0 * std::cos(config.omega * time);

            return {Ex, 0.0, 0.0};
        }
    }
    return {};
}


Vec3 magnetic_field(
    const Vec3& position,
    double time,
    const FieldConfig& config
)
{
    (void)position;
    (void)time;

    switch (config.type) {

        case FieldCase::UniformB:
            // B = B0 z-hat
            return {0.0, 0.0, config.B0};


        case FieldCase::ExB:
            // Uniform magnetic field:
            // B = B0 z-hat
            return {0.0, 0.0, config.B0};


        case FieldCase::GradB: {
            // Weakly nonuniform magnetic field:
            //
            // B(x) = B0 * (1 + gradient * x) z-hat
            //
            // This field satisfies div(B) = 0 because
            // Bx = By = 0 and Bz has no z dependence.
            const double Bz = config.B0 * (1.0 + config.gradient * position.x);
            return {0.0, 0.0, Bz};
        }


        case FieldCase::Curvature: {
            // Constant-magnitude azimuthal magnetic field:
            //
            // B = B0 * e_phi
            //
            // In cartesian coordinates:
            //
            // e_phi = (-y/R, x/R, 0)
            // where R = sqrt(x^2 + y^2)
            //
            // The field magnitude is constant, so there is no gradB drift
            // The curved field lnes produce curvature drift
            const Vec3 radial_position{
                position.x,
                position.y,
                0.0
            };
            const double radius = radial_position.norm();

            if (radius < 1.0e-12) {
                throw std::runtime_error("Curvature field is singular at R = 0.");
            };

            return {
                -config.B0 * position.y / radius,
                config.B0 * position.x / radius,
                0.0
            };
        }


        case FieldCase::NonuniformE:
            // Uniform magnetic field:
            //
            // B = B0 z-hat

            return {0.0, 0.0, config.B0};


        case FieldCase::Polarization:
            // Uniform magnetic field:
            //
            // B = B0 z-hat

            return {0.0, 0.0, config.B0};
    }

    return {};
}
