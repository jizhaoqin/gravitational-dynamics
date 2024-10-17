#include "../../const-SI.h"
#include "../../utils.h"
#include <cmath>
#include <iostream>
#include <math.h>

#define TMAX (1000 * YEAR / TUNIT)
#define SunMass (MSUN / MUNIT)
#define EarthMass (MEARTH / MUNIT)
#define EarthA (1. * AU / LUNIT)

const double dt = 0.1 * YEAR / TUNIT;

void initializer(state &Target, double Scaler = 1.2) {
  Target.time = 0;
  Target.mass = EarthMass;
  Target.pos = Vector3(1. * EarthA, 0, 0);
  Target.vel = Vector3(0, Scaler * 2 * M_PI * EarthA, 0);
}

double totalEnergy(state &Target, double MassStar) {
  double T = 0.5 * MassStar * Target.mass * Target.vel.modulusSquare();
  double U = -G * MassStar * Target.mass / Target.pos.modulus();
  return T + U;
}

void updateAcc(state &Target) {
  double distance = Target.pos.modulus();
  Target.acc = (-G * SunMass / pow(distance, 3)) * Target.pos;
}

void updateVel(state &Target, double dt) { Target.vel = Target.vel + Target.acc * dt; }

void updatePos(state &Target, double dt) { Target.pos = Target.pos + Target.vel * dt; }

void output(state &Target) { std::cout << Target.time << ' ' << Target.pos << Target.vel << totalEnergy(Target, SunMass) << std::endl; }

//
//
//
//
////////************************* BEGIN main function *************************////////
int main(int argc, const char **argv) {
  state EarthStatus;
  initializer(EarthStatus, 1.2);

  int OrbitCounter = 0;
  while (EarthStatus.time < TMAX) {
    output(EarthStatus);

    // Euler method
    // updateAcc(EarthStatus);
    // updateVel(EarthStatus, dt);
    // updatePos(EarthStatus, dt);

    // DKD
    updatePos(EarthStatus, 0.5 * dt);
    updateAcc(EarthStatus);
    updateVel(EarthStatus, 1.0 * dt);
    updatePos(EarthStatus, 0.5 * dt);

    EarthStatus.time += dt;
  }
}
////////************************* END of main function *************************////////
