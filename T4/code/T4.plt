
  set term postscript enhanced color solid rounded
  set output 'T4.ps'
  set encoding iso_8859_1
  set colors classic

# Note: make sure the labels correspond to the data we're plotting!
  set title 'start velocity = 0.8 v_{circ}'

  iEul = 1; tEul = 'Euler, stepsize 0.0025'
  iMid = 2; tMid = 'Midpoint rule, stepsize 0.005'
  iRK4 = 3; tRK4 = 'RK 4, stepsize 0.01'
  iSym = 4; tSym = 'symplectic Euler, stepsize 0.0025'
  iDKD = 5; tDKD = 'DKD Leapfrog, stepsize 0.0025'
  iKDK = 6; tKDK = 'KDK Leapfrog, stepsize 0.0025'
  iSun = 7

  set style line iEul lt 1 pt 7 ps .2 lc rgb '#ff0000'
  set style line iMid lt 1 pt 7 ps .2 lc rgb '#00cc00'
  set style line iRK4 lt 1 pt 7 ps .2 lc rgb '#0000ff'
  set style line iSym lt 1 pt 7 ps .2 lc rgb '#ffdd00' lw 1.5
  set style line iKDK lt 1 pt 7 ps .2 lc rgb '#ff00ff'
  set style line iDKD lt 1 pt 7 ps .2 lc rgb '#00ffff' lw .5
  set style line iSun lt 1 pt 7 ps 1. lc rgb '#000000'

  o(i,j,k,x) = (int(i)%j==k)?(x):(1/0)

  set grid lt 1 lc rgb '#eeeeee'
  set tics front

#-----------------------------------------------------------------------
# Data files are assumed to have the following columns:
#
#     time  xpos ypos zpos  xvel yvel zvel  Ekin Epot  orbit#
#
# Nr:   1     2    3    4     5    6    7     8    9     10
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# plot #1: orbit

  r = 5
  set xrange [-r:r]
  set yrange [-r:r]
  set size ratio -1
  set xlabel 'x'
  set ylabel 'y'
  set key top left reverse Left opaque invert

  plot \
    'T4.sym.out' u (o($10,100,0,$2)):3 w l ls iSym t tSym.' (every 100th orbit)', \
    'T4.kdk.out' u (o($10,100,0,$2)):3 w l ls iKDK t tKDK.' (every 100th orbit)', \
    'T4.dkd.out' u (o($10,100,0,$2)):3 w l ls iDKD t tDKD.' (every 100th orbit)', \
    'T4.rk4.out' u (o($10,100,0,$2)):3 w l ls iRK4 t tRK4.' (every 100th orbit)', \
    'T4.mid.out' u (o($10, 50,0,$2)):3 w l ls iMid t tMid.' (every 50th orbit)', \
    'T4.eul.out' u (o($10,  1,0,$2)):3 w l ls iEul t tEul.' (every orbit)', \
    '-' u 1:2 w p ls iSun notitle
    0 0
    e

#-----------------------------------------------------------------------
# plot #2: orbit zoom

  r = 1.5
  set xrange [-r:r]
  set yrange [-r:r]
  replot

#-----------------------------------------------------------------------
# plot #3: energy vs. time, linear time

  set auto x
  set auto y
  set size noratio
  set xlabel 't'
  set ylabel 'E/m'
  set key top right noreverse Right

  plot \
    'T4.sym.out' u 1:($8+$9) w l ls iSym t tSym, \
    'T4.kdk.out' u 1:($8+$9) w l ls iKDK t tKDK, \
    'T4.dkd.out' u 1:($8+$9) w l ls iDKD t tDKD, \
    'T4.rk4.out' u 1:($8+$9) w l ls iRK4 t tRK4, \
    'T4.mid.out' u 1:($8+$9) w l ls iMid t tMid, \
    'T4.eul.out' u 1:($8+$9) w l ls iEul t tEul

#-----------------------------------------------------------------------
# plot #4: energy vs. time, logarithmic time

  set log x
  set key top left reverse Left
  replot

#-----------------------------------------------------------------------

