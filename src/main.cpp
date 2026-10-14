#include <iostream>

#include "Vec3.hpp"
#include "Particle.hpp"
#include "Fields.hpp"
#include "Integrator.hpp"


int main()
{
    std::cout << "Charged Particle Simulator\n";


    // -------------------------
    // Test 1: cross product
    // -------------------------

    Vec3 ex{1.0, 0.0, 0.0};
    Vec3 ey{0.0, 1.0, 0.0};

    Vec3 ez = cross(ex, ey);

    std::cout
        << "ex x ey = "
        << ez.x << " "
        << ez.y << " "
        << ez.z << '\n';


    // -------------------------
    // Test 2: particle
    // -------------------------

    Particle particle{
        {0.0, 0.0, 0.0},   // position
        {1.0, 0.0, 1.0},   // velocity
        1.0,                // mass
        1.0                 // charge
    };


    // -------------------------
    // Test 3: uniform B
    // -------------------------

    FieldConfig fields;

    fields.type = FieldCase::UniformB;
    fields.B0 = 1.0;


    Vec3 B = magnetic_field(
        particle.position,
        0.0,
        fields
    );

    std::cout
        << "B = "
        << B.x << " "
        << B.y << " "
        << B.z << '\n';


    // -------------------------
    // Test 4: acceleration
    // -------------------------

    Vec3 a = acceleration(
        particle,
        particle.position,
        particle.velocity,
        0.0,
        fields
    );

    std::cout
        << "a = "
        << a.x << " "
        << a.y << " "
        << a.z << '\n';


    return 0;
}
