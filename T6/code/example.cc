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

#define NTEST 50
#define MASS (msun / m_unit)

using namespace std;

#include "../../const.h"
#include "example.h"

// For variable time steps:
const double dt_factor_time_dependent = 0.1;

double RAD_PLUMMER = 0;

void initialize_sol_system(vector<state> &start) {
  start.push_back(state(0.3871 * au / l_unit, 47.87e5 / v_unit, 3.302e26 / m_unit)); // Mercury
  start.push_back(state(0.723 * au / l_unit, 35.02e5 / v_unit, 4.869e27 / m_unit));  // Venus
  start.push_back(state(1.0 * au / l_unit, 29.78e5 / v_unit, 5.97219e27 / m_unit));  // Earth
  start.push_back(state(1.524 * au / l_unit, 24.13e5 / v_unit, 6.219e26 / m_unit));  // Mars
  start.push_back(state(5.203 * au / l_unit, 13.07e5 / v_unit, 1.899e30 / m_unit));  // Jupiter
  start.push_back(state(9.5826 * au / l_unit, 9.69e5 / v_unit, 5.685e29 / m_unit));  // Saturn
  start.push_back(state(19.201 * au / l_unit, 6.81e5 / v_unit, 8.683e29 / m_unit));  // Uranus
  start.push_back(state(30.047 * au / l_unit, 5.43e5 / v_unit, 1.0243e29 / m_unit)); // Neptune
  start.push_back(state(39.482 * au / l_unit, 4.72e5 / v_unit, 1.25e25 / m_unit));   // Pluto
}

void initialize_2(vector<state> &start) {
  start.push_back(state(au / l_unit, 0.02 * 29.78e5 / v_unit, 1));
  start.push_back(state(au / l_unit, 0.1 * 29.78e5 / v_unit, 1));
}

// Set up a "chain" of test particles.
void initialize_chain(vector<state> &start) {
  for (int i = 0; i < NTEST; i++) {
    state test(au / l_unit, 29.78e5 / v_unit * (float)i / NTEST, 1);
    test.x.comp[1] = 0.0 + (float)i / NTEST / 10 * au / l_unit;
    start.push_back(test);
  }
}

double effective_potential(state &here) {
  vec L = here.x.cross(here.v);
  double pot = G * here.m * MASS / sqrt(here.x.abs2() + RAD_PLUMMER * RAD_PLUMMER); /* Potential energy */
  double red = 0.5 * here.m * L.abs2() / here.x.abs2();
  return (red - pot);
}

double total_energy(state &here) {
  double pot = here.m * MASS * G / here.x.abs(); /* Potential energy */
  double kin = 0.5 * here.m * here.v.abs2();     /* Kinetic energy */
  return (kin - pot);                            /* Here we return the total energy */
}

void calculate_acceleration(vector<state> &here) {
  for (int i = 0; i < here.size(); i++) {
    double r = sqrt(here[i].x.abs2() + RAD_PLUMMER * RAD_PLUMMER);
    here[i].a = here[i].x * (-G * MASS / r / r / r);
  }
}

// Update velocities
void kick(vector<state> &here, double mydt) {
  for (int i = 0; i < here.size(); i++)
    here[i].v += here[i].a * mydt;
}

// Update positions
void drift(vector<state> &here, double mydt) {
  for (int i = 0; i < here.size(); i++)
    here[i].x += here[i].v * mydt;
}

////////************************* BEGIN main function *************************////////
int main(int argc, const char **argv) {
  /* State is a structure (defined it in example.h) which contains all quantities (like position, velocity, etc.)
     We create a vector of it which we call "current" and update it for each time step. */

  vector<state> current;       // Vector of particles
  double next_output_time = 0; // Time for next output
  double dt = 0;               // Timestep to be defined later
  ofstream outfile;            // File to directly output the data to

  // Check if an argument was passed to our program.
  if (argc < 2) {
    cout << "please give a number to the program:" << endl
         << "  1: Solar system" << endl
         << "  2: Two particles" << endl
         << "  3: Particle chain" << endl;
    return (1);
  }

  // Initialize the system depending on the user's choice.
  switch (atoi(argv[1])) {
  case 1:
    initialize_sol_system(current);
    // RAD_PLUMMER = 0; // for point mass.
    RAD_PLUMMER = (0.5 * au / l_unit);
    outfile.open("sol_system.dat");
    break;
  case 2:
    initialize_2(current);
    RAD_PLUMMER = (0.5 * au / l_unit);
    outfile.open("plummer_2.dat");
    break;
  case 3:
    initialize_chain(current);
    RAD_PLUMMER = (0.5 * au / l_unit);
    outfile.open("plummer_chain.dat");
    break;
  default:
    cout << "Only 1, 2, or 3 is acceptable!" << endl;
    return (2);
  }

  // useful information
  outfile << "# L = " << l_unit << endl;
  outfile << "# T = " << t_unit << endl;
  outfile << "# M = " << m_unit << endl;
  outfile << "# V = " << v_unit << endl;
  outfile << "# G = " << G << endl;

  // initial acceleration
  calculate_acceleration(current);

  ////////************************* LOOP BEGIN  *************************////////
  while (current[0].t < TMAX) {
    // output every OUTPUT_FREQUENCY
    if (current[0].t >= next_output_time) {
      cout << current[0].t << " " << dt << " " << next_output_time << endl;

      outfile << current[0].t << " ";
      for (int i = 0; i < current.size(); i++)
        outfile << current[i].x << current[i].v << " " << effective_potential(current[i]) << " " << total_energy(current[i]) << " ";

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
////////************************* END of main function *************************////////

/* How to get it running:
  1) make
  2) ./myprogram
*/
