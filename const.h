// cgs units

// Some useful constants
#define year 3.1557600e7    // year in s
#define msun 1.9884e33      // Sun mass in g
#define mmoon 7.3477e25     // Moon mass in g
#define mearth 5.9722e27    // Earth mass in g
#define mjupiter 1.898e30   // Jupiter mass in g
#define au 1.495978707e13   // astronomical unit in cm
#define rearth 6371e5       // Earth radius in cm
#define dmoon 384400e5      // Earth-Moon distance in cm
#define parsec 3.085678e18  // parsec in cm
#define kparsec 3.085678e21 // kiloparsec in cm
#define mparsec 3.085678e24 // megaparsec in cm

// unit system of the simulation (solar system)
#define m_unit (msun)
#define l_unit (au)
#define t_unit (year)

// unit system of the simulation (galaxies)
// #define m_unit (1e10*msun)
// #define l_unit (kparsec)
// #define t_unit (1e6*year) // megayears

// unit system of the simulation (clusters)
// #define m_unit (1e10*msun)
// #define l_unit (kparsec)
// #define t_unit (1e9*year) // gigayears

// resulting velocity unit
#define v_unit (l_unit / t_unit)

// Gravitational constant and speed of light in internal units
#define G (6.6743e-8 * m_unit * t_unit * t_unit / l_unit / l_unit / l_unit)
#define c (2.99792458e+10 * t_unit / l_unit)
