#include "Fields.hpp"


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
            return {};


        case FieldCase::NonuniformE:
            return {};


        case FieldCase::Polarization:
            return {};
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


        case FieldCase::Curvature:
            return {};


        case FieldCase::NonuniformE:
            return {};


        case FieldCase::Polarization:
            return {};
    }

    return {};
}
