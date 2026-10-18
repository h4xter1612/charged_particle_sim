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

    const double q_over_m =
        particle.charge / particle.mass;

    const Vec3 lorentz_term =
        E + cross(velocity, B);

    return q_over_m * lorentz_term;
}

void rk4_step(
    Particle& particle,
    const FieldConfig& fields,
    double time,
    double dt
)
{
    const Vec3 r0 = particle.position;
    const Vec3 v0 = particle.velocity;

    // k1: derivatives at the beginning of the time step.
    const Vec3 k1_r = v0;
    const Vec3 k1_v = acceleration(
        particle,
        r0,
        v0,
        time,
        fields
    );

    // k2: derivatives at the midpoint estimated from k1.
    const Vec3 r2 =
        r0 + k1_r * (0.5 * dt);
    const Vec3 v2 =
        v0 + k1_v * (0.5 * dt);

    const Vec3 k2_r = v2;
    const Vec3 k2_v = acceleration(
        particle,
        r2,
        v2,
        time + 0.5 * dt,
        fields
    );

    // k3: derivatives at the midpoint estimated from k2.
    const Vec3 r3 =
        r0 + k2_r * (0.5 * dt);
    const Vec3 v3 =
        v0 + k2_v * (0.5 * dt);

    const Vec3 k3_r = v3;
    const Vec3 k3_v = acceleration(
        particle,
        r3,
        v3,
        time + 0.5 * dt,
        fields
    );

    // k4: derivatives at the end of the time step estimated from k3.
    const Vec3 r4 =
        r0 + k3_r * dt;
    const Vec3 v4 =
        v0 + k3_v * dt;

    const Vec3 k4_r = v4;
    const Vec3 k4_v = acceleration(
        particle,
        r4,
        v4,
        time + dt,
        fields
    );

    // Fourth-order weighted combination.
    particle.position =
        r0
        + (dt / 6.0)
        * (
            k1_r
            + 2.0 * k2_r
            + 2.0 * k3_r
            + k4_r
        );

    particle.velocity =
        v0
        + (dt / 6.0)
        * (
            k1_v
            + 2.0 * k2_v
            + 2.0 * k3_v
            + k4_v
        );
}
