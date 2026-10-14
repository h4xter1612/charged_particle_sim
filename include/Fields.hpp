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

    // Lo usaremos despues para gradientes
    double gradient{0.0};

    // Lo usaremos para campos dependientes del tiempo
    double omega{0.0};
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
