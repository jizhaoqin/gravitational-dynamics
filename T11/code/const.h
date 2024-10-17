
// Some useful constants
#define year 3.1557600e7    // year in s
#define msun 1.9884e33      // Sun mass in g
#define mmoon 7.342e25      // Moon mass in g
#define mearth 5.9722e27    // Earth mass in g
#define au 1.495978707e13   // astronomical unit in cm
#define rearth 6371e5       // Earth radius in cm
#define dmoon 384400e5      // Earth-Moon distance in cm
#define parsec 3.085678e18  // parsec in cm
#define kparsec 3.085678e21 // kiloparsec in cm
#define mparsec 3.085678e24 // megaparsec in cm

// unit system of the simulation (solar system)
// #define m_unit (msun)
// #define l_unit (au)
// #define t_unit (year)
// #define v_unit (l_unit / t_unit) // resulting velocity unit

// unit system of the simulation (galaxies)
#define m_unit (1e10 * msun)
#define l_unit (kparsec)
#define v_unit (1e5)             // 1e5 cm/s = 1km/s
#define t_unit (l_unit / v_unit) // resulting time unit

// unit system of the simulation (clusters)
// #define m_unit (1e10*msun)
// #define l_unit (kparsec)
// #define t_unit (1e9*year)        // gigayears
// #define v_unit (l_unit / t_unit) // resulting velocity unit

// Gravitational constant and speed of light in internal units
#define G (6.6743e-8 * m_unit * t_unit * t_unit / l_unit / l_unit / l_unit)
#define c (2.99792458e+10 * t_unit / l_unit)
