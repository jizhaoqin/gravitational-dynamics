#include <cmath>
#include <iostream>

// We can set TMAX either here or in the Makefile.
#ifndef TMAX
#define TMAX 20
#endif

using namespace std;

#include "const.h"
#include "example.h"

// Time integration parameter:
const double dt = 0.01;

// Set up the earth orbit (1 au, 1 year period) for the state which we start with.
void initialize(state &start) {
  start.t = 0;
  start.m = mearth / m_unit;
  start.x = vec(au / l_unit, 0, 0);
  start.v = vec(0, 2.0 * M_PI, 0);
  // Eccentric orbit
  // start.v = vec(0,0.8 * 2.0 * M_PI,0);
}

/* Calculate total energy at a given state. This is a function which gets a state as an argument
   and returns a double (the total energy of the state). */
double calculate_total_energy(state &here) {
  double mass = msun / m_unit;                   /* Mass of sun in internal units
                                                    (will be 1 in the chosen units) */
  double pot = here.m * mass * G / here.x.abs(); /* Potential energy */
  double kin = 0.5 * here.m * here.v.abs2();     /* Kinetic energy */

  return (kin - pot); /* Here we return the total energy */
}

// Calculate acceleration of a given state. This is a procedure which gets a state as an argument.
void calculate_acceleration(state &here) {
  double r = here.x.abs();
  double sol_mass = msun / m_unit; /* Mass of sun in internal units
                                      (will be 1 in the chosen units) */
  here.a = here.x * (-G * sol_mass / r / r / r);
}

// Update velocities of a given state. This is a procedure which gets a state as an argument.
void update_velocity(state &here) { here.v += here.a * dt; }

// Update positions of a given state. This is a procedure which gets a state as an argument.
void update_position(state &here) { here.x += here.v * dt; }

int main(int argc, const char **argv) {
  /* State is a structure (we defined it in example.h) which contains all quantities (like position, velocity, etc.)
     needed to describe a test particle.
     We create one instance of it which we call "current" and update it for each time step. */

  state current;

  // Output at the beginning some useful information ...
  cout << "# L = " << l_unit << endl;
  cout << "# T = " << t_unit << endl;
  cout << "# M = " << m_unit << endl;
  cout << "# V = " << v_unit << endl;
  cout << "# G = " << G << endl;

  // Initialize the system. We defined this procedure before.
  initialize(current);

  // This is loop which is repeated until integration of t reaches our defined TMAX.
  while (current.t < TMAX) {
    // Output time, position, and total energy of the current time.
    cout << current.t << " " << current.x << " " << calculate_total_energy(current) << endl;

    // Now we calculate the acceleration of the current state.
    calculate_acceleration(current);

    /* We first update the position of the current state.
       Remember, here the velocity is used, so if we would
       first update the velocity this would be different! */
    update_position(current);

    // Now we update the velocity of the current state.
    update_velocity(current);

    // Finally we update the time t.
    current.t += dt;
  }

  return (0);
}

/* How to get it running:
  1) make
  2) ./myprogram > orbit.dat
  3) gnuplot
     gnuplot> plot 'orbit.dat' u ($2):($3)
*/
