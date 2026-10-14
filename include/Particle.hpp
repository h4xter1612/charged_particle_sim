#pragma once

#include "Vec3.hpp"

struct Particle {
    Vec3 position;
    Vec3 velocity;

    double mass;
    double charge;

    double kinetic_energy() const {
        // TODO:
        // K = 1/2 m v^2
        //
        // Ya tendras velocity.dot(velocity)
        // o velocity.norm()

        return 0.5 * mass * velocity.dot(velocity);
    }
};
