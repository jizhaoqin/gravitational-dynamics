#include <math.h>
#include <ostream>

/* Vector class for working with vectors in three-dimensional space. */
class vec {
public:
  double comp[3]; // The components of the vector.

  // Define some standard ways of initializing,
  // which allows us to write:
  //   vec x(1,0,0);

  vec() { comp[0] = comp[1] = comp[2] = 0; }
  vec(double _x, double _y, double _z) {
    comp[0] = _x;
    comp[1] = _y;
    comp[2] = _z;
  }

  // Define some standard functions one would like to have for vectors,
  // which allows to do:
  //   double len=x.abs();

  double abs() { return (sqrt(comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2])); }
  double abs2() { return (comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2]); }
  double dot(const vec &_x) { return (comp[0] * _x.comp[0] + comp[1] * _x.comp[1] + comp[2] * _x.comp[2]); }

  // Cross product can be written as:
  //   vec z=x.cross(y);
  vec cross(const vec &_x) {
    return (vec(comp[1] * _x.comp[2] - comp[2] * _x.comp[1], comp[2] * _x.comp[0] - comp[0] * _x.comp[2],
                comp[0] * _x.comp[1] - comp[1] * _x.comp[0]));
  }

  // Define some simple mathematical operators,
  // which allows us to write things like:
  //   vec y=x*3;
  vec operator*(double _x) { return (vec(comp[0] * _x, comp[1] * _x, comp[2] * _x)); }
  vec operator/(double _x) { return (vec(comp[0] / _x, comp[1] / _x, comp[2] / _x)); }
  vec operator+(const vec &_x) { return (vec(comp[0] + _x.comp[0], comp[1] + _x.comp[1], comp[2] + _x.comp[2])); }
  vec operator-(const vec &_x) { return (vec(comp[0] - _x.comp[0], comp[1] - _x.comp[1], comp[2] - _x.comp[2])); }

  vec &operator+=(const vec &_x) {
    this->comp[0] += _x.comp[0];
    this->comp[1] += _x.comp[1];
    this->comp[2] += _x.comp[2];
    return *this;
  }
  vec &operator-=(const vec &_x) {
    this->comp[0] -= _x.comp[0];
    this->comp[1] -= _x.comp[1];
    this->comp[2] -= _x.comp[2];
    return *this;
  }
};

inline vec operator*(double _a, const vec &_x) { return (vec(_a * _x.comp[0], _a * _x.comp[1], _a * _x.comp[2])); }

// Define the format for printing vectors:

inline std::ostream &operator<<(std::ostream &os, const vec &_x) {
  os << ' ' << _x.comp[0] << ' ' << _x.comp[1] << ' ' << _x.comp[2] << ' ';
  return os;
}

/* Here we define the variables of the system which we
   want to keep for each time within one structure. */

struct state {
  vec x, v, a;
  double m, t;

  state() {
    x = v = a = vec(0, 0, 0);
    m = t = 0;
  }
  state(double _x, double _v, double _m) {
    x = v = a = vec(0, 0, 0);
    t = 0;
    x.comp[0] = _x;
    v.comp[1] = _v;
    m = _m;
  }
  state(double _x, double _y, double _z, double _vx, double _vy, double _vz, double _m) {
    x = vec(_x, _y, _z);
    v = vec(_vx, _vy, _vz);
    t = 0;
    m = _m;
  }
};
