#include <omp.h>
#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// We can set TMAX either here or in the Makefile.
#ifndef TMAX
#define TMAX 10
#endif

#ifndef OUTPUT_FREQUENCY
#define OUTPUT_FREQUENCY 0.01
#endif

using namespace std;

#include "const.h"
#include "example.h"

// For variable time steps:
#ifdef USE_TIMESTEP_LEVELS
const double dt_factor_time_dependent = 0.5;
const double dt_base = 0.000001;
#else
const double dt_factor_time_dependent = 2.5;
const double dt_base = 1;
#endif

double SOFT = 0.1;

// Read initial conditions from file, no change
void initialize_read_xvm(vector<state> &start, string name, double x0, double v0)
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

            start.push_back(state(x + x0, y, z, vx + v0, vy, vz, m));
        }
        infile.close();
    }
    else
        cout << "Unable to open file: " << name << endl;
    cout << "# Read " << start.size() << " points !" << endl;
}

// Read initial conditions from file, exchange x and z coordinates
void initialize_read_vertical(vector<state> &start, string name, double x0, double v0)
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
            // start.push_back(state(x + x0, y, z, vx + v0, vy, vz, m));
            start.push_back(state(z, y, x + x0, vx + v0, vy, vz, m)); // exchange x and z
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
void calculate_acceleration(vector<state> &here, vector<int> myactive)
{
#pragma omp parallel for
    for (int j = 0; j < myactive.size(); j++)
    {
        int i = myactive[j];
        here[i].a = vec();
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
void kick_active(vector<state> &here, double mydt, vector<int> myactive)
{
#pragma omp parallel for
    for (int j = 0; j < myactive.size(); j++)
    {
        int i = myactive[j];
        here[i].v += here[i].a * here[i].dt * mydt;
    }
}

// Update velocities of a given state. This is a procedure which gets a vector of state and a deltaT as arguments.
void kick_all(vector<state> &here, double mydt)
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

// Re-calculate the timesteps for the active particles.
void set_timesteps(vector<state> &here, vector<int> myactive)
{
#pragma omp parallel for
    for (int j = 0; j < myactive.size(); j++)
    {
        int i = myactive[j];
        here[i].t += here[i].dt;
        here[i].dt = pow(2, floor(log(dt_factor_time_dependent / here[i].a.abs() / dt_base) / log(2)));

        if (here[i].dt < 1)
        {
            cout << "particle " << i << " reached too small timestep " << dt_factor_time_dependent / here[i].a.abs() << ", decrease timebase " << dt_base << endl;
            std::terminate();
        }
    }
}

// Create the list of active particles.
int find_next_active(vector<state> &here, MY_TIME_TYPE mytime, vector<int> &myactive)
{
    int dt_next = std::numeric_limits<int>::max();

    for (int i = 0; i < here.size(); i++)
    {
        int dt_this = here[i].t + here[i].dt - mytime;
        if (dt_this < dt_next && dt_this > 0)
            dt_next = here[i].t + here[i].dt - mytime;
    }
    myactive.resize(0);
    for (int i = 0; i < here.size(); i++)
        if (here[i].t + here[i].dt - mytime == dt_next)
            myactive.push_back(i);

    if (dt_next <= 0)
    {
        cout << "dt_next " << dt_next << " for " << myactive.size() << " particles, should not happen !!! " << endl;
        std::terminate();
    }

    return dt_next;
}

// Set the next timestep.
double find_min_timestep(vector<state> &here)
{
    double a_max = 0;
    for (int i = 0; i < here.size(); i++)
        if (here[i].a.abs() > a_max)
            a_max = here[i].a.abs();

    return dt_factor_time_dependent / a_max;
}

// #################################  main function  ################################# //
int main(int argc, const char **argv)
{
    int split_files = 0;                    // Flag to output individual files
    vector<int> active;                     // List of active particles
    vector<state> current;                  // Vector of our particle system
    double next_output_time = 0;            // Time for next output
    MY_TIME_TYPE dt = 0, time = 0;          // System time
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
             << "  6: read initial conditions from file" << endl
             << "  7: read initial conditions from two files and displace them by dx and dvx" << endl;
        return (1);
    }

    // Initialize the system depending on the user's choice.
    switch (atoi(argv[1]))
    {
    case 6:
        if (argc < 3)
            initialize_read_xvm(current, "plummer_1000.ascii", 0, 0);
        else
            initialize_read_xvm(current, argv[2], 0, 0);
        add_massive_particle(current);
        split_files = 0;
        break;
    case 7:
        if (argc < 3)
        {
            initialize_read_vertical(current, "./data/v160_c12_l05_1e3.ic.ascii", 0, 0);
            initialize_read_xvm(current, "./data/v160_c12_l05_1e3.ic.ascii", 50, -140);
        }
        else
        {
            if (argc < 6)
            {
                cout << "Please use: ./myprogram 7 <file1> <file2> dx dvx" << endl;
                return (0);
            }
            else
            {
                initialize_read_xvm(current, argv[2], 0, 0);
                initialize_read_xvm(current, argv[3], atof(argv[4]), atof(argv[5]));
            }
        }
        split_files = 0;
        if (split_files == 0)
            outfile.open("./data/merge-1e3-vertical.dat");
        break;
    default:
        cout << "Only 1, 2, 3, 4, 5, 6, or 7 are acceptable!" << endl;
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

    // Now we calculate the acceleration of the initial state, all particles are set to active.
    for (int i = 0; i < current.size(); i++)
        active.push_back(i);
    calculate_acceleration(current, active);
#ifdef USE_TIMESTEP_LEVELS
    set_timesteps(current, active);
    dt = find_next_active(current, time, active);
#else
    dt = find_min_timestep(current);
#endif

    // #################################  main loop  ################################# //
    while (time * dt_base < TMAX)
    {
        if (time * dt_base >= next_output_time)
        {
            if (split_files == 0)
            {
                cout << time << " " << next_output_time << endl;
                outfile << time * dt_base << " ";
                for (int i = 0; i < current.size(); i++)
                {
                    outfile << current[i].x << " ";
                }
                outfile << endl;
            }
            else
            {
                char name[100];
                sprintf(name, "pos_%04d.dat", count_outputs);
                cout << "Writing " << name << " t = " << time * dt_base << " Nstep = " << count_steps << endl;
                outfile.open(name);
                outfile << "# " << time * dt_base << " " << count_outputs << endl;
                for (int i = 0; i < current.size(); i++)
                    outfile << current[i].x << endl;
                outfile.close();
            }
            count_outputs++;
            next_output_time += OUTPUT_FREQUENCY;
        }

        // KDK Scheme:
#ifdef USE_TIMESTEP_LEVELS
        kick_active(current, 0.5 * dt_base, active);  // Kick only active particle
        drift(current, dt * dt_base);                 // Drift system to next time
        calculate_acceleration(current, active);      // Acceleration only for active particle
        kick_active(current, 0.5 * dt_base, active);  // Kick only active particle
        time += dt;                                   // Advance system time
        set_timesteps(current, active);               // Calculate new timestep of active particle
        dt = find_next_active(current, time, active); // find next list of active particles
#else
        kick_all(current, 0.5 * dt);
        drift(current, dt);                      // Drift system to next time
        calculate_acceleration(current, active); // Acceleration only for active particle
        kick_all(current, 0.5 * dt);
        time += dt;                      // Advance system time
        dt = find_min_timestep(current); // find minimal timestep
#endif
        count_steps++;
    }

    if (split_files == 0)
    {
        outfile.close();
    }

    return (0);
}

/* How to get it running:
  1) make
  2) ./myprogram 7
*/
