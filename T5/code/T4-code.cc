#include "../../const.h"
#include "../../example.h"
// #include <cmath>
#include <iostream>

#ifndef TMAX
#define TMAX 250
#endif

using namespace std;

const double dt_factor = 0.005;
const double dt_factor_time_dependent = 0.5;
const int OUTPUT_FREQUENCY = 40;

void initialize(state &start);
double totalEnergy(state &here);
void updateAcc(state &here);
void updateVel(state &here, double mydt);
void updatePos(state &here, double mydt);

int main(int argc, const char **argv) {
  state current;

  // Output at the beginning some useful information ...
  cout << "# L = " << l_unit << endl;
  cout << "# T = " << t_unit << endl;
  cout << "# M = " << m_unit << endl;
  cout << "# V = " << v_unit << endl;
  cout << "# G = " << G << endl;

  int orbit_counter = 0;
  double dt = 20 * dt_factor; // here the unit is YEAR/TUNIT = 1

  initialize(current);
  updateAcc(current);

  while (current.t < TMAX) {
    // Output time, position, and total energy of the current time.
    // if (orbit_counter % OUTPUT_FREQUENCY == 0) {
    //   cout << current.t << current.x << current.v << totalEnergy(current) << endl;
    // }
    cout << current.t << current.x << current.v << totalEnergy(current) << endl;

    // for orbit counter
    // double y_old = current.x.comp[1];

    // Variable timestep:
    // dt = dt_factor_time_dependent / current.a.abs();

    // KDK Scheme:
    // updateVel(current, 0.5 * dt);
    // updatePos(current, 1.0 * dt);
    // updateAcc(current);
    // updateVel(current, 0.5 * dt);

    // DKD Scheme:
    updatePos(current, 0.5 * dt);
    updateAcc(current);
    updateVel(current, 1.0 * dt);
    updatePos(current, 0.5 * dt);

    // to count the orbits (as sign reversal from - to + in the y position):
    // if (y_old < 0 && current.x.comp[1] >= 0)
    //   orbit_counter++;

    current.t += dt;
  }

  return (0);
}

void initialize(state &start) {
  start.t = 0;
  start.m = mearth / m_unit;
  start.x = vec(au / l_unit, 0, 0);
  start.v = vec(0, 1.2 * 2.0 * M_PI, 0); // elliptical or circular
}

double totalEnergy(state &here) {
  double mass = msun / m_unit;
  double pot = here.m * mass * G / here.x.abs();
  double kin = 0.5 * here.m * here.v.abs2();
  return (kin - pot);
}

void updateAcc(state &here) {
  double r = here.x.abs();
  double sol_mass = msun / m_unit;
  here.a = here.x * (-G * sol_mass / r / r / r);
}

void updateVel(state &here, double mydt) { here.v += here.a * mydt; }

void updatePos(state &here, double mydt) { here.x += here.v * mydt; }