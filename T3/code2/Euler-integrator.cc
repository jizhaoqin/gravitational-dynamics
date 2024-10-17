#include "../../const-SI.h"
#include "../../utils.h"
#include <iostream>

#define SunMass (MSUN / MUNIT)
#define EarthMass (MEARTH / MUNIT)
#define EarthA (1. * AU / LUNIT)

void test();
void updateStatus(double DT, Vector3 &Position, Vector3 &Velocity, Vector3 &Acceleration);
double getEnergy(double Mass1, double Mass2, double Distance, Vector3 Velocity);
void output(int i, double DT, double MassPlanet, double MassStar, Vector3 Position, Vector3 Velocity);

int main(int argc, const char **argv) {
  /**
   * initial conditions of the Earth
   * a is dependent on position, so only two independent status,
   * and it would make different which one is updated first.
   */
  double circularVelocity = sqrt(G * SunMass / EarthA);
  Vector3 EarthPosition(EarthA, 0, 0);
  Vector3 EarthVelocity(0, 1.2 * circularVelocity, 0);
  Vector3 EarthAcceleration(0, 0, 0);

  // update status and print out data
  int Step = 250;
  double DT = 1e-1 * YEAR / TUNIT;
  for (int i = 0; i < Step; i++) {
    updateStatus(DT, EarthPosition, EarthVelocity, EarthAcceleration);
    output(i, DT, EarthMass, SunMass, EarthPosition, EarthVelocity);
  }
}

// update position first or velocity first would lead to totally different results
void updateStatus(double DT, Vector3 &Position, Vector3 &Velocity, Vector3 &Acceleration) {
  double Distance = Position.modulus();                      // calculate distance
  Vector3 RHat = Position.normalize() * (-1);                  // calculate direction vector
  Acceleration = RHat * (G * SunMass / (Distance * Distance)); // calculate acceleration

  Velocity = Velocity + Acceleration * DT; // update velocity
  Position = Position + Velocity * DT;     // update position
}

double getEnergy(double Mass1, double Mass2, double Distance, Vector3 Velocity) {
  double T = 0.5 * Mass1 * Velocity.modulus() * Velocity.modulus();
  double U = -G * Mass2 * Mass1 / Distance;
  return T + U;
}

void output(int i, double DT, double MassPlanet, double MassStar, Vector3 Position, Vector3 Velocity) {
  double Distance = Position.modulus();
  std::cout << (i + 1) * DT << ' ';
  std::cout << Position.x << ' ' << Position.y << ' ' << Position.z << ' ';
  std::cout << getEnergy(MassPlanet, MassStar, Distance, Velocity) << std::endl;
}

// check unit is desired
void test() {
  std::cout << SunMass << '\t' << EarthA << '\t' << G << std::endl;
  std::cout << sqrt(G * SunMass / EarthA) << std::endl;
  std::cout << EarthMass << std::endl;
}