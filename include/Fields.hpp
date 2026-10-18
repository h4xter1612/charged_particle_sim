#pragma once

#include "Vec3.hpp"

enum class FieldCase {
    UniformB,
    ExB,
    GradB,
    Curvature,
    NonuniformE,
    Polarization
};

struct FieldConfig {
    FieldCase type{FieldCase::UniformB};

    double B0{1.0};
    double E0{0.0};

    // Spatial coefficient in Bz(x) = B0 * (1 + gradient * x).
    double gradient{0.0};

    // Angular frequency of the time-dependent electric field.
    double omega{0.0};

    // Spatial wavenumber of the nonuniform electric field.
    double wavenumber{0.0};
};

Vec3 electric_field(
    const Vec3& position,
    double time,
    const FieldConfig& config
);

Vec3 magnetic_field(
    const Vec3& position,
    double time,
    const FieldConfig& config
);
