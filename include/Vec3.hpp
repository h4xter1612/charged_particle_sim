#pragma once

#include <cmath>
#include <stdexcept>

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

    Vec3 operator-(const Vec3& other) const {
        return {
            x - other.x,
            y - other.y,
            z - other.z
        };
    }

    Vec3 operator*(double scalar) const {
        return {
            x * scalar,
            y * scalar,
            z * scalar
        };
    }

    Vec3 operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::invalid_argument("Cannot divide Vec3 by zero.");
        }

        return {
            x / scalar,
            y / scalar,
            z / scalar
        };
    }

    double dot(const Vec3& other) const {
        return
            x * other.x
            + y * other.y
            + z * other.z;
    }

    double norm() const {
        return std::sqrt(dot(*this));
    }
};

inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

inline Vec3 operator*(double scalar, const Vec3& v)
{
    return v * scalar;
}
