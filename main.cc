#include "const.h"
#include "subroutine.h"
#include "utils.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>

// We can set TMAX either here or in the Makefile.
#ifndef TMAX
#define TMAX 2000
#endif

#ifndef OUTPUT_FREQUENCY
#define OUTPUT_FREQUENCY 0.1
#endif

using namespace std;

// For variable time steps:
const double dt_factor_time_dependent = 0.001;

double SOFT = 0;

// double distance2SmallBody(double R, double M0, double M1);
// void initialize_sem(vector<state> &start);
// void initialize_L123(vector<state> &start, double pert_x, double pert_y);
// void initialize_L4(vector<state> &start, double pert_x, double pert_y);
// void initialize_swing(vector<state> &start);
// void calculate_acceleration_fixed(vector<state> &here);
// void calculate_acceleration(vector<state> &here);
// void kick(vector<state> &here, double mydt);
// void drift(vector<state> &here, double mydt);

////////************************* BEGIN main function*************************////////
int main(int argc, const char **argv) {
    vector<state> current;       // particle system,target
    double next_output_time = 0; // Time for next output
    double dt = 0;               // Timestep
    ofstream outfile;

    // Check if an argument was passed to our program.
    if (argc < 2) {
        cout << "please give a number to the program:" << endl
             << "  1: Test particle in L123, give perturbation in x and y as extra arguments" << endl
             << "  2: Sun-Earth-Moon system" << endl
             << "  3: Swing-by" << endl
             << "  4: Test particle in L45, give perturbation in x and y as extra arguments!" << endl;
        return (1);
    }

    // Initialize the system depending on the user's choice.
    switch (atoi(argv[1])) {
    case 1:
        if (argc < 3)
            initialize_L123(current, 0, 0);
        else if (argc < 4)
            initialize_L123(current, atof(argv[2]), 0);
        else
            initialize_L123(current, atof(argv[2]), atof(argv[3]));
        outfile.open("./data/L123.dat");
        break;
    case 2:
        initialize_sem(current);
        outfile.open("./data/earth_moon_sun.dat");
        break;
    case 3:
        initialize_swing(current);
        outfile.open("./data/swing.dat");
        break;
    case 4:
        if (argc < 3)
            initialize_L4(current, 0, 0);
        else if (argc < 4)
            initialize_L4(current, atof(argv[2]), 0);
        else
            initialize_L4(current, atof(argv[2]), atof(argv[3]));
        outfile.open("./data/L4.dat");
        break;
    default:
        cout << "Only 1, 2, 3, or 4 are acceptable!" << endl;
        return (2);
    }

    // Output at the beginning some useful information ...
    outfile.precision(10);
    outfile << "# L = " << l_unit << endl;
    outfile << "# T = " << t_unit << endl;
    outfile << "# M = " << m_unit << endl;
    outfile << "# V = " << v_unit << endl;
    outfile << "# G = " << G << endl;

    // Now we calculate the acceleration of the initial state.
    calculate_acceleration(current);

    // ************************* LOOP BEGIN  ************************* //
    while (current[0].t < TMAX) {
        // Write output (trick: do it only if time defined by OUTPUT_FREQUENCY
        // is reached!).
        if (current[0].t >= next_output_time) {
            cout << current[0].t << " " << dt << " " << next_output_time << endl;

            outfile << current[0].t << " ";
            for (int i = 0; i < current.size(); i++)
                outfile << current[i].x << current[i].v << " ";
            outfile << endl;
            next_output_time += OUTPUT_FREQUENCY;
        }

        // Variable timestep (trick: first find largest acceleration and
        // convert to timestep at the end!).
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
////////************************* END main function*************************////////
