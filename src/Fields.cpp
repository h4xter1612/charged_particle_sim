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
            return {
                config.E0,
                0.0,
                0.0
            };


        case FieldCase::GradB:
            // Para este caso probablemente E = 0
            return {};


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
            return {
                0.0,
                0.0,
                config.B0
            };


        case FieldCase::GradB:
            // TODO mas adelante.
            //
            // Aqui tendremos MUCHO cuidado con:
            //
            // div B = 0
            //
            return {};


        case FieldCase::Curvature:
            return {};


        case FieldCase::NonuniformE:
            return {};


        case FieldCase::Polarization:
            return {};
    }

    return {};
}
