#ifndef CONSTANTS
#define CONSTANTS

// Some useful constants
#define YEAR 3.1557600e7    // YEAR in s
#define MSUN 1.9884e30      // Sun mass in kg
#define MMOON 7.342e22      // Moon mass in kg
#define MEARTH 5.9722e24    // Earth mass in kg
#define AU 1.495978707e11   // astronomical unit in m
#define REARTH 6371e3       // Earth radius in m
#define DMOON 384400e3      // Earth-Moon distance in m
#define PARSEC 3.085678e16  // PARSEC in m
#define KPARSEC 3.085678e19 // kiloPARSEC in m
#define MPARSEC 3.085678e22 // megaPARSEC in m

// unit system of the simulation (solar system)
#define MUNIT (MSUN)
#define LUNIT (AU)
#define TUNIT (YEAR)

// unit system of the simulation (galaxies)
// #define MUNIT (1e10*MSUN)
// #define LUNIT (kPARSEC)
// #define TUNIT (1e6*YEAR) // megayears

// unit system of the simulation (clusters)
// #define MUNIT (1e10*MSUN)
// #define LUNIT (kPARSEC)
// #define TUNIT (1e9*YEAR) // gigayears

// resulting velocity unit
#define vUNIT (LUNIT / TUNIT)
#define aUNIT (LUNIT / TUNIT / TUNIT)

// Gravitational constant and speed of light in internal units
#define G (6.6743e-11 * MUNIT * TUNIT * TUNIT / LUNIT / LUNIT / LUNIT)
#define c (2.99792458e+8 * TUNIT / LUNIT)

#endif