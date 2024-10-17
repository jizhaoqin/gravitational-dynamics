#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>

// We can set TMAX either here or in the Makefile.
#ifndef TMAX
#define TMAX 20
#endif

#ifndef OUTPUT_FREQUENCY
#define OUTPUT_FREQUENCY 0.01
#endif

using namespace std;

#include "../../const.h"
#include "example.h"

// For variable time steps:
const double dt_factor_time_dependent = 0.1;

double SOFT = 0;

// Set up the solar system: we set up the Sun and several planets which we collect into a "vector".
void initialize_sol_system(vector<state> &start) {
  start.push_back(state(0 / l_unit, 0 / v_unit, msun / m_unit));                     // Sun
  start.push_back(state(0.3871 * au / l_unit, 47.87e5 / v_unit, 3.302e26 / m_unit)); // Mercury
  start.push_back(state(0.723 * au / l_unit, 35.02e5 / v_unit, 4.869e27 / m_unit));  // Venus
  start.push_back(state(1.0 * au / l_unit, 29.78e5 / v_unit, mearth / m_unit));      // Earth
  start.push_back(state(1.524 * au / l_unit, 24.13e5 / v_unit, 6.219e26 / m_unit));  // Mars
  start.push_back(state(5.203 * au / l_unit, 13.07e5 / v_unit, 1.899e30 / m_unit));  // Jupiter
  start.push_back(state(9.5826 * au / l_unit, 9.69e5 / v_unit, 5.685e29 / m_unit));  // Saturn
  start.push_back(state(19.201 * au / l_unit, 6.81e5 / v_unit, 8.683e29 / m_unit));  // Uranus
  start.push_back(state(30.047 * au / l_unit, 5.43e5 / v_unit, 1.0243e29 / m_unit)); // Neptune
  start.push_back(state(39.482 * au / l_unit, 4.72e5 / v_unit, 1.25e25 / m_unit));   // Pluto
}

// Set up the Sun-Earth-Moon system.
void initialize_ems(vector<state> &start) {
  double m0 = msun / m_unit;
  double m1 = mearth / m_unit;
  double m2 = 7.349e22 / m_unit;

  start.push_back(state(-m1 / (m0 + m1) * au / l_unit, -29.78e5 * m1 / m0 / v_unit, m0));
  start.push_back(state(m0 / (m0 + m1) * au / l_unit, 29.78e5 / v_unit, m1));
  start.push_back(state(start[1].x.comp[0] - 384.4e8 / l_unit, 29.78e5 / v_unit + 1.023e5 / v_unit, m2));
}

// Set up swing-by.
void initialize_swing(vector<state> &start) {
  double m0 = msun / m_unit;
  double m1 = mearth / m_unit * 10000; // Make Earth more massive
  double m2 = 1e6 / m_unit;            // Just a small value

  start.push_back(state(-m1 / (m0 + m1) * au / l_unit, -29.78e5 * m1 / m0 / v_unit, m0));
  start.push_back(state(m0 / (m0 + m1) * au / l_unit, 29.78e5 / v_unit, m1));

  double d_phi = 18;
  double r_0 = 0.1;
  double v_0 = 129.0e5 / v_unit;

  state test;
  test.x.comp[0] = r_0 * cos(d_phi / 180 * M_PI);
  test.x.comp[1] = r_0 * sin(d_phi / 180 * M_PI);
  test.v.comp[0] = v_0 * cos(d_phi / 180 * M_PI);
  test.v.comp[1] = v_0 * sin(d_phi / 180 * M_PI);
  test.m = m2;
  start.push_back(test);

  test.x.comp[0] = r_0 * cos((d_phi - 1) / 180 * M_PI);
  test.x.comp[1] = r_0 * sin((d_phi - 1) / 180 * M_PI);
  test.v.comp[0] = v_0 * cos((d_phi - 1) / 180 * M_PI);
  test.v.comp[1] = v_0 * sin((d_phi - 1) / 180 * M_PI);
  start.push_back(test);

  test.x.comp[0] = r_0 * cos((d_phi + 1) / 180 * M_PI);
  test.x.comp[1] = r_0 * sin((d_phi + 1) / 180 * M_PI);
  test.v.comp[0] = v_0 * cos((d_phi + 1) / 180 * M_PI);
  test.v.comp[1] = v_0 * sin((d_phi + 1) / 180 * M_PI);
  start.push_back(test);
}

// Calculate acceleration of a given state. This is a procedure which gets a vector of state as an argument.
void calculate_acceleration(vector<state> &here) {
  for (int i = 0; i < here.size(); i++) {
    here[i].a = vec(); // zero all components
    for (int j = 0; j < here.size(); j++)
      if (i != j) {
        vec dist_vec = here[i].x - here[j].x;
        double r = sqrt(dist_vec.abs2() + SOFT);
        here[i].a += -G * here[j].m / (r * r * r) * dist_vec;
      }
  }
}

// Update velocities of a given state. This is a procedure which gets a vector of state and a deltaT as arguments.
void kick(vector<state> &here, double mydt) {
  for (int i = 0; i < here.size(); i++)
    here[i].v += here[i].a * mydt;
}

// Update positions of a given state. This is a procedure which gets a vector of state and a deltaT as arguments.
void drift(vector<state> &here, double mydt) {
  for (int i = 0; i < here.size(); i++)
    here[i].x += here[i].v * mydt;
}

////////************************* BEGIN main function *************************////////
int main(int argc, const char **argv) {
  /* State is a structure (we defined it in example.h) which contains all quantities (like position, velocity, etc.)
     needed to describe a test particle.
     We create a vector of it which we call "current" and update it for each time step. */

  vector<state> current;       // Vector of our particle system
  double next_output_time = 0; // Time for next output
  double dt = 0;               // Timestep to be defined later
  ofstream outfile;            // File to directly output the data to

  // Check if an argument was passed to our program.
  if (argc < 2) {
    cout << "please give a number to the program:" << endl
         << "  1: Solar system" << endl
         << "  2: Sun-Earth-Moon system" << endl
         << "  3: Swing-by" << endl;
    return (1);
  }

  // Initialize the system depending on the user's choice.
  switch (atoi(argv[1])) {
  case 1:
    initialize_sol_system(current);
    outfile.open("sol_system.dat");
    break;
  case 2:
    initialize_ems(current);
    outfile.open("earth_moon_sun.dat");
    break;
  case 3:
    initialize_swing(current);
    outfile.open("swing.dat");
    break;
  default:
    cout << "Only 1, 2, or 3 is acceptable!" << endl;
    return (2);
  }

  // Output at the beginning some useful information ...
  outfile << "# L = " << l_unit << endl;
  outfile << "# T = " << t_unit << endl;
  outfile << "# M = " << m_unit << endl;
  outfile << "# V = " << v_unit << endl;
  outfile << "# G = " << G << endl;

  // Now we calculate the acceleration of the initial state.
  calculate_acceleration(current);

  // This is loop which is repeated until integration of t reaches our defined TMAX.
  while (current[0].t < TMAX) {
    // Write output (trick: do it only if time defined by OUTPUT_FREQUENCY is reached!).
    if (current[0].t >= next_output_time) {
      cout << current[0].t << " " << dt << " " << next_output_time << endl;

      outfile << current[0].t << " ";
      for (int i = 0; i < current.size(); i++)
        outfile << current[i].x << current[i].v << " ";
      outfile << endl;
      next_output_time += OUTPUT_FREQUENCY;
    }

    // Variable timestep (trick: first find largest acceleration and convert to timestep at the end!).
    double a_max = 0;
    for (int i = 0; i < current.size(); i++)
      if (current[i].a.abs() > a_max)
        a_max = current[i].a.abs();
    dt = dt_factor_time_dependent / a_max;

    // KDK Scheme:
    kick(current, 0.5 * dt);
    drift(current, 1.0 * dt);
    calculate_acceleration(current);
    kick(current, 0.5 * dt);

    // Finally we update the time t.
    for (int i = 0; i < current.size(); i++)
      current[i].t += dt;
  }

  outfile.close();

  return (0);
}
////////************************* END main function *************************////////

/* How to get it running:
  1) make
  2) ./myprogram
*/
