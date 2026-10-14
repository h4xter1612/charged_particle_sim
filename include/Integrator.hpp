#pragma once

#include "Particle.hpp"
#include "Fields.hpp"


Vec3 acceleration(
    const Particle& particle,
    const Vec3& position,
    const Vec3& velocity,
    double time,
    const FieldConfig& fields
);


void rk4_step(
    Particle& particle,
    const FieldConfig& fields,
    double time,
    double dt
);


// Lo dejamos declarado desde ahora,
// pero NO lo implementamos todavia.
void boris_step(
    Particle& particle,
    const FieldConfig& fields,
    double time,
    double dt
);
