#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdlib.h>
#include <vector>

// We can set TMAX either here or in the Makefile.
#ifndef TMAX
#define TMAX 200
#endif

#ifndef OUTPUT_FREQUENCY
#define OUTPUT_FREQUENCY 0.005
#endif

using namespace std;

#include "const.h"
#include "example.h"

// For variable time steps:
const double dt_factor_time_dependent = 30;

double SOFT = 0.01;

void initialize_sol_system(vector<state> &start)
{
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

  // Change reference frame so that total momentum is zero so that the system does not drift!
  double m = 0;
  vec mom(0, 0, 0);
  for (int i = 0; i < start.size(); i++)
  {
    m += start[i].m;
    mom += start[i].v * start[i].m;
  }
  mom = mom / m;
  for (int i = 0; i < start.size(); i++)
    start[i].v -= mom;
}

void initialize_ems(vector<state> &start)
{
  double m0 = msun / m_unit;
  double m1 = mearth / m_unit;
  double m2 = 7.349e22 / m_unit;

  start.push_back(state(-m1 / (m0 + m1) * au / l_unit, -29.78e5 * m1 / m0 / v_unit, m0));
  start.push_back(state(m0 / (m0 + m1) * au / l_unit, 29.78e5 / v_unit, m1));
  start.push_back(state(start[1].x.comp[0] - 384.4e8 / l_unit, 29.78e5 / v_unit + 1.023e5 / v_unit, m2));
}

// Set up planetary swing-by.
void initialize_swing(vector<state> &start)
{
  double m0 = msun / m_unit;
  double m1 = mearth / m_unit * 10000; // Make Earth more massive
  double m2 = 1e6 / m_unit;            // Just a small value

  start.push_back(state(-m1 / (m0 + m1) * au / l_unit, -29.78e5 * m1 / m0 / v_unit, m0));
  start.push_back(state(m0 / (m0 + m1) * au / l_unit, 29.78e5 / v_unit, m1));

  double d_phi = 34;
  double r_0 = 0.1;
  double v_0 = 129.0e5 / v_unit;
  double phase[3] = {0, -10, 2.5};

  state test;
  test.m = m2;
  for (int i = 0; i < 3; i++)
  {
    test.x.comp[0] = r_0 * cos((d_phi + phase[i]) / 180 * M_PI);
    test.x.comp[1] = r_0 * sin((d_phi + phase[i]) / 180 * M_PI);
    test.v.comp[0] = v_0 * cos((d_phi + phase[i]) / 180 * M_PI);
    test.v.comp[1] = v_0 * sin((d_phi + phase[i]) / 180 * M_PI);
    start.push_back(test);
  }
}

// Set up L4.
void initialize_l4(vector<state> &start, double pert_x, double pert_y)
{
  double m0 = msun / m_unit;
  // double m1 = msun / m_unit;
  double m1 = mearth / m_unit;
  double m2 = 1e8 / m_unit; // test body

  // Set up the sun-planet system with zero momentum to keep it from drifting away!
  double r = au / l_unit;                                        // distance between planet and sun
  double v = sqrt(G * m0 / (r * (1 + m1 / m0)));                 // orbital velocity of planet
  start.push_back(state(-m1 / (m0 + m1) * r, -v * m1 / m0, m0)); // reduced velocity of sun
  start.push_back(state(m0 / (m0 + m1) * r, v, m1));             // velocity of planet

  // The L4 point is always defined by an equilateral triangle with the sun-planet system.
  state test;
  test.m = m2;
  test.x.comp[0] = start[0].x.comp[0] + cos(60. / 180 * M_PI) * r;
  test.x.comp[1] = sin(60. / 180 * M_PI) * r;
  // Velocity of L4 is given by its distance from the center of mass
  // and the orbital period of the sun-planet system.
  double vt = test.x.abs() / start[1].x.abs() * v;
  test.v.comp[0] = -test.x.comp[1] / test.x.abs() * vt;
  test.v.comp[1] = test.x.comp[0] / test.x.abs() * vt;

  // Add some perturbation in case.
  test.x.comp[0] += pert_x;
  test.x.comp[1] += pert_y;

  start.push_back(test);
}

// Set up a homogeneous sphere of particles.
void initialize_sphere(vector<state> &start)
{
#define NGRID 5
  double m = 1.0;

  for (int ix = -NGRID; ix <= NGRID; ix++)
    for (int iy = -NGRID; iy <= NGRID; iy++)
      for (int iz = -NGRID; iz <= NGRID; iz++)
      {
        vec pos(ix, iy, iz);
        state test;
        test.x = pos * au / l_unit;
        test.m = m;
        if (test.x.abs() <= NGRID * au / l_unit)
          start.push_back(test);
      }
  cout << "# Created " << start.size() << " points !" << endl;
}

// Read initial conditions from file.
void initialize_read_xvm(vector<state> &start, string name)
{
  string line;
  ifstream infile;

  infile.open(name.c_str());

  if (infile.is_open())
  {
    while (getline(infile, line))
    {
      stringstream ss(line);
      double x, y, z, vx, vy, vz, m;
      ss >> x >> y >> z >> vx >> vy >> vz >> m;

      start.push_back(state(x, y, z, vx, vy, vz, m));
    }
    infile.close();
  }
  else
  {
    cout << "Unable to open file: " << name << endl;
  }
  cout << "# Read " << start.size() << " points !" << endl;
}

// Add extra particle.
void add_massive_particle(vector<state> &start)
{
  state test;
  test.m = start[0].m * 20; // Test Particle with 20 times more mass
  test.x.comp[0] = 5.0;     // Place Test Particle at r=5
  double mtot = 0;
  for (int i = 0; i < start.size(); i++)
    if (start[i].x.abs() < test.x.comp[0])
      mtot += start[i].m;

  test.v.comp[1] = sqrt(G * mtot / test.x.comp[0]);
  cout << "Orbit and Orbital velocity of test particle = " << test.x.comp[0] << "," << test.v.comp[1] << endl;
  cout << "G*M/vc and log(Lambda): " << G * test.m / test.v.comp[1] << "," << log(mtot / test.m) << endl;
  cout << "Factor: " << 2 * 0.428 * log(mtot / test.m) * G * test.m / test.v.comp[1] << endl;

  start.push_back(test);

  // Change reference frame so that total momentum is zero so that the system does not drift!
  mtot = 0;
  vec mom(0, 0, 0);
  vec pos(0, 0, 0);
  for (int i = 0; i < start.size(); i++)
  {
    mtot += start[i].m;
    mom += start[i].v * start[i].m;
    pos += start[i].x * start[i].m;
  }
  mom = mom / mtot;
  pos = pos / mtot;
  for (int i = 0; i < start.size(); i++)
  {
    start[i].v -= mom;
    start[i].x -= pos;
  }
}

// Calculate acceleration of a given state. This is a procedure which gets a vector of state as an argument.
void calculate_acceleration_fixed(vector<state> &here)
{
  for (int i = 1; i < here.size(); i++)
  {
    vec dist_vec = here[i].x - here[0].x;
    double r = sqrt(dist_vec.abs2() + SOFT);
    here[i].a = -G * here[0].m / (r * r * r) * dist_vec;
  }
}

// Calculate acceleration of a given state. This is a procedure which gets a vector of state as an argument.
void calculate_acceleration(vector<state> &here)
{
#pragma omp parallel for
  for (int i = 0; i < here.size(); i++)
  {
    here[i].a = vec(); // zero all components
    for (int j = 0; j < here.size(); j++)
      if (i != j)
      {
        vec dist_vec = here[i].x - here[j].x;
        double r = sqrt(dist_vec.abs2() + SOFT);
        here[i].a += -G * here[j].m / (r * r * r) * dist_vec;
      }
  }
}

// Update velocities of a given state. This is a procedure which gets a vector of state and a deltaT as arguments.
void kick(vector<state> &here, double mydt)
{
#pragma omp parallel for
  for (int i = 0; i < here.size(); i++)
    here[i].v += here[i].a * mydt;
}

// Update positions of a given state. This is a procedure which gets a vector of state and a deltaT as arguments.
void drift(vector<state> &here, double mydt)
{
#pragma omp parallel for
  for (int i = 0; i < here.size(); i++)
    here[i].x += here[i].v * mydt;
}

////////************************* Begin main function *************************////////
int main(int argc, const char **argv)
{

  int split_files = 0;                    // Flag to output individual files
  vector<state> current;                  // Vector of our particle system
  double next_output_time = 0;            // Time for next output
  double dt = 0;                          // Timestep to be defined later
  ofstream outfile;                       // File to directly output the data to
  int count_outputs = 0, count_steps = 0; // Counter for output files

#pragma omp parallel
  {
#pragma omp master
    {
      cout << "Using " << omp_get_num_threads() << " OpenMP threads ..." << endl;
    }
  }

  // Check if an argument was passed to our program.
  if (argc < 2)
  {
    cout << "please give a number to the program:" << endl
         << "  1: Solar system" << endl
         << "  2: Sun-Earth-Moon system" << endl
         << "  3: Swing-by" << endl
         << "  4: Test particle in L4, give perturbation in x and y as extra arguments!" << endl
         << "  5: homogeneous density sphere for collapse test" << endl
         << "  6: read initial conditions from file" << endl;
    return (1);
  }

  // Initialize the system depending on the user's choice.
  switch (atoi(argv[1]))
  {
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
  case 4:
    if (argc < 3)
      initialize_l4(current, 0, 0);
    else if (argc < 4)
      initialize_l4(current, atof(argv[2]), 0);
    else
      initialize_l4(current, atof(argv[2]), atof(argv[3]));
    outfile.open("l4.dat");
    break;
  case 5:
    initialize_sphere(current);
    outfile.open("sphere.dat");
    break;
  case 6:
    if (argc < 3)
      initialize_read_xvm(current, "./data/plummer_1000.ascii");
    else
      initialize_read_xvm(current, argv[2]);
    add_massive_particle(current);

    split_files = 0; // 0: into one file, 1: into mutiple files
    if (split_files == 0)
      outfile.open("./data/friction.dat");
    break;
  default:
    cout << "Only 1, 2, 3, 4, 5, or 6 are acceptable!" << endl;
    return (2);
  }

  // Output at the beginning some useful information ...
  if (split_files == 0)
  {
    outfile << "# L = " << l_unit << endl;
    outfile << "# T = " << t_unit << endl;
    outfile << "# M = " << m_unit << endl;
    outfile << "# V = " << v_unit << endl;
    outfile << "# G = " << G << endl;
    outfile << "# N = " << current.size() << endl;
  }
  else
  {
    cout << "# L = " << l_unit << endl;
    cout << "# T = " << t_unit << endl;
    cout << "# M = " << m_unit << endl;
    cout << "# V = " << v_unit << endl;
    cout << "# G = " << G << endl;
    cout << "# N = " << current.size() << endl;
  }

  // Now we calculate the acceleration of the initial state.
  calculate_acceleration(current);

  // ************************* LOOP BEGIN  ************************* //
  while (current[0].t < TMAX)
  {
    // Write output (trick: do it only if time defined by OUTPUT_FREQUENCY is reached!).
    if (current[0].t >= next_output_time)
    {
      if (split_files == 0)
      { // append to one file
        cout << current[0].t << " " << dt << " " << next_output_time << endl;

        outfile << current[0].t << " ";
        for (int i = 0; i < current.size(); i++)
          outfile << current[i].x << " ";
        outfile << endl;
      }
      else
      { // make individual files
        char name[100];
        sprintf(name, "./data/pos_%04d.dat", count_outputs);
        cout << "Writing " << name << " t = " << current[0].t << ", dt= " << dt << " Nstep = " << count_steps << endl;
        outfile.open(name);
        outfile << "# " << current[0].t << " " << count_outputs << endl;
        for (int i = 0; i < current.size(); i++)
          outfile << current[i].x << endl;
        outfile.close();
      }
      count_outputs++;
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

    count_steps++;
  }

  if (split_files == 0)
    outfile.close();

  return (0);
}

/* How to get it running:
  1) make
  2) ./myprogram 6
*/
