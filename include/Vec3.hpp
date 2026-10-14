#pragma once

#include <cmath>

struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};

    Vec3 operator+(const Vec3& other) const {
        return {
            x + other.x,
            y + other.y,
            z + other.z
        };
    }

    // TODO 1:
    // Implementa resta vectorial:
    // a - b
    Vec3 operator-(const Vec3& other) const {
        // TODO
        return {
            x - other.x,
            y - other.y,
            z - other.z
        };
    }

    // TODO 2:
    // Multiplicacion por escalar:
    // v * scalar
    Vec3 operator*(double scalar) const {
        // TODO
        return {
            x * scalar,
            y * scalar,
            z * scalar
        };
    }

    // TODO 3:
    // Division por escalar:
    Vec3 operator/(double scalar) const {
        // TODO
        return {
            x / scalar,
            y / scalar,
            z / scalar

        };
    }

    // TODO 4:
    // Producto punto
    double dot(const Vec3& other) const {
        // TODO
        return (x * other.x) + (y * other.y) + (z * other.z);
    }

    // TODO 5:
    // Norma |v|
    double norm() const {
        // Pista:
        // puedes usar dot(*this)
        return std::sqrt(dot(*this));
    }
};


// TODO 6:
// Producto cruz
inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    // Recuerda:
    //
    // a x b =
    // ( ay*bz - az*by,
    //   az*bx - ax*bz,
    //   ax*by - ay*bx )

    return {
        (a.y * b.z) - (a.z * b.y),
        (a.z * b.x) - (a.x * b.z),
        (a.x * b.y) - (a.y * b.x)
    };
}


// Esta te la doy porque luego resulta comodo escribir:
//
// 2.0 * v
//
// ademas de:
//
// v * 2.0
//
inline Vec3 operator*(double scalar, const Vec3& v)
{
    return v * scalar;
}
