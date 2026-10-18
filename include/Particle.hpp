#pragma once

#include "Vec3.hpp"

struct Particle {
    Vec3 position;
    Vec3 velocity;

    double mass;
    double charge;

    double kinetic_energy() const {
        return 0.5 * mass * velocity.dot(velocity);
    }
};
