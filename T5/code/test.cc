#include "../../const-SI.h"
#include "../../utils.h"
#include <cmath>
#include <iostream>
#include <math.h>

void testPrint(const Vector3 &Vector);
void testConstants();

int main(int argc, const char **argv) {
  Vector3 DemoVector(1., 2., 3.);
  // testPrint(DemoVector);
  testConstants();

  return 0;
}

void testPrint(const Vector3 &Vector) {
  std::cout << Vector << std::endl;
  std::cout << Vector.dot(Vector) << std::endl;
  std::cout << (Vector * 3) << std::endl;
}

void testConstants() {
  std::cout << G << std::endl;
  std::cout << sqrt(G * MSUN / MUNIT / pow((AU / LUNIT), 3)) << std::endl;
}
