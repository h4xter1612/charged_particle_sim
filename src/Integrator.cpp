#include "Integrator.hpp"


Vec3 acceleration(
    const Particle& particle,
    const Vec3& position,
    const Vec3& velocity,
    double time,
    const FieldConfig& fields
)
{
    const Vec3 E = electric_field(
        position,
        time,
        fields
    );

    const Vec3 B = magnetic_field(
        position,
        time,
        fields
    );

    // TODO:
    //
    // Implementa:
    //
    // a = (q/m) * (E + v x B)
    //
    // Ya tendras disponibles:
    //
    // cross(...)
    // operator+
    // operator*
    //
    Vec3 a;
    const double q_over_m = particle.charge / particle.mass;
    const Vec3 lorentz_term = E + cross(velocity, B);

    return q_over_m * lorentz_term;
}


void rk4_step(
    Particle& particle,
    const FieldConfig& fields,
    double time,
    double dt
)
{
    // TODO:
    //
    // Todavia NO lo hagas.
    //
    // Primero vamos a derivar juntos
    // RK4 para el sistema:
    //
    // dr/dt = v
    // dv/dt = a(r,v,t)
}


void boris_step(
    Particle& particle,
    const FieldConfig& fields,
    double time,
    double dt
)
{
    // TODO futuro

    (void)particle;
    (void)fields;
    (void)time;
    (void)dt;
}
