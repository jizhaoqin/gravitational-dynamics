#ifndef MYUTILS_HPP
#define MYUTILS_HPP

#include <cmath>
#include <cstdio>
#include <iostream>

class Vector3 {
public:
  double x, y, z;

  // Default constructor
  Vector3() : x(0), y(0), z(0) {}

  // Constructor with initial values
  Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  // Copy constructor
  Vector3(const Vector3 &other) : x(other.x), y(other.y), z(other.z) {}

  // Print the vector
  void print() const { printf("(%.4f, %.4f, %.4f)\n", x, y, z); }

  // Assignment operator
  Vector3 &operator=(const Vector3 &other) {
    if (this != &other) {
      x = other.x;
      y = other.y;
      z = other.z;
    }
    return *this;
  }

  // Addition operator
  Vector3 operator+(const Vector3 &other) const { return Vector3(x + other.x, y + other.y, z + other.z); }

  // Subtraction operator
  Vector3 operator-(const Vector3 &other) const { return Vector3(x - other.x, y - other.y, z - other.z); }

  // Multiplication by scalar operator
  Vector3 operator*(double scalar) const { return Vector3(x * scalar, y * scalar, z * scalar); }

  // Division by scalar operator
  Vector3 operator/(double scalar) const { return Vector3(x / scalar, y / scalar, z / scalar); }

  // Dot product
  double dot(const Vector3 &other) const { return x * other.x + y * other.y + z * other.z; }

  // Cross product
  Vector3 cross(const Vector3 &other) const {
    return Vector3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
  }

  // Magnitude
  double modulus() const { return sqrt(x * x + y * y + z * z); }
  double modulusSquare() const { return x * x + y * y + z * z; }

  // Normalize
  Vector3 normalize() const {
    double mag = modulus();
    if (mag > 0) {
      return *this / mag;
    } else {
      return Vector3();
    }
  }
};

// Overload the * operator inline to allow scalar multiplication commutativity
inline Vector3 operator*(double scaler, const Vector3 &v) { return (Vector3(scaler * v.x, scaler * v.y, scaler * v.z)); }

// Overload the << operator to print vectors
inline std::ostream &operator<<(std::ostream &os, const Vector3 &v) {
  os << v.x << ' ' << v.y << ' ' << v.z << ' ';
  return os;
}

// status of a single particle
struct state {
  Vector3 pos, vel, acc;
  double mass, time;
};

#endif
