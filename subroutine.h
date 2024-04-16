#include "const.h"
#include "utils.h"
#include <vector>

using namespace std;

// R is the distance between the 2 major bodies
inline double distance2SmallBody(double R, double M0, double M1) { return R * pow(M1 / 3 / M0, 1. / 3.); }

// Sun-Earth-Moon system
inline void initialize_sem(vector<state> &start) {
    double m0 = msun / m_unit;
    double m1 = mearth / m_unit;
    double m2 = 7.349e22 / m_unit;

    start.push_back(state(-m1 / (m0 + m1) * au / l_unit, -29.78e5 * m1 / m0 / v_unit, m0));
    start.push_back(state(m0 / (m0 + m1) * au / l_unit, 29.78e5 / v_unit, m1));
    start.push_back(state(start[1].x.comp[0] - 384.4e8 / l_unit, 29.78e5 / v_unit + 1.023e5 / v_unit, m2));
}

// Jupiter at Earth's position, L123 points setup
inline void initialize_L123(vector<state> &start, double pert_x, double pert_y) {
    double m0 = msun / m_unit;
    // double m1 = mjupiter / m_unit;
    double m1 = mearth / m_unit;
    double m2 = 1e6 / m_unit; // test body

    // Set up the sun-planet system with zero momentum to keep it from drifting away!
    double r = au / l_unit;                                        // distance between planet and sun
    double v = sqrt(G * m0 / (r * (1 + m1 / m0)));                 // orbital velocity of planet
    start.push_back(state(-m1 / (m0 + m1) * r, -v * m1 / m0, m0)); // reduced velocity of sun
    start.push_back(state(m0 / (m0 + m1) * r, v, m1));             // velocity of planet

    // take the assumption: m1 << m0
    state test;
    test.m = m2;
    // double distance = distance2SmallBody(r, m0, m1); // L12
    // test.x.comp[0] = start[1].x.comp[0] - distance; // L1
    // test.x.comp[0] = start[1].x.comp[0] + distance; // L2
    // double v_total = test.x.abs() / start[1].x.abs() * v;
    // test.v.comp[1] = v_total;

    double distance = r * 7 * m1 / 12 / m0;                // L3
    test.x.comp[0] = start[1].x.comp[0] * (-1) + distance; // L3
    double v_total = test.x.abs() / start[1].x.abs() * v;  // L3
    test.v.comp[1] = v_total * (-1);                       // L3

    // Add some perturbation in case.
    test.x.comp[0] += pert_x;
    test.x.comp[1] += pert_y;

    start.push_back(test);
}

// Set up planetary swing-by.
inline void initialize_swing(vector<state> &start) {
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
    for (int i = 0; i < 3; i++) {
        test.x.comp[0] = r_0 * cos((d_phi + phase[i]) / 180 * M_PI);
        test.x.comp[1] = r_0 * sin((d_phi + phase[i]) / 180 * M_PI);
        test.v.comp[0] = v_0 * cos((d_phi + phase[i]) / 180 * M_PI);
        test.v.comp[1] = v_0 * sin((d_phi + phase[i]) / 180 * M_PI);
        start.push_back(test);
    }
}

// Set up L4.
inline void initialize_L4(vector<state> &start, double pert_x, double pert_y) {
    double m0 = msun / m_unit;
    double m1 = mjupiter / m_unit;
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
    // total velocity proportional to the distance to COM, components could be derived by the similarity of triangles
    double v_total = test.x.abs() / start[1].x.abs() * v;
    test.v.comp[0] = -test.x.comp[1] / test.x.abs() * v_total;
    test.v.comp[1] = test.x.comp[0] / test.x.abs() * v_total;

    // Add some perturbation in case.
    test.x.comp[0] += pert_x;
    test.x.comp[1] += pert_y;

    start.push_back(test);
}

// with a single central force field
inline void calculate_acceleration_fixed(vector<state> &here, double SOFT = 0) {
    for (int i = 1; i < here.size(); i++) {
        vec dist_vec = here[i].x - here[0].x;
        double r = sqrt(dist_vec.abs2() + SOFT);
        here[i].a = -G * here[0].m / (r * r * r) * dist_vec;
    }
}

// interactions between each two particles
inline void calculate_acceleration(vector<state> &here, double SOFT = 0) {
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
inline void kick(vector<state> &here, double mydt) {
    for (int i = 0; i < here.size(); i++)
        here[i].v += here[i].a * mydt;
}

// Update positions of a given state. This is a procedure which gets a vector of state and a deltaT as arguments.
inline void drift(vector<state> &here, double mydt) {
    for (int i = 0; i < here.size(); i++)
        here[i].x += here[i].v * mydt;
}