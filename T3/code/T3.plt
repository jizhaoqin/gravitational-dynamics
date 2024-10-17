
  set term postscript enhanced color solid
  set output 'T3.ps'
  set encoding iso_8859_1
  set colors classic

  set grid lt 1 lc rgb '#eeeeee'
  set tics front

#-----------------------------------------------------------------------
# The data file "T3.out" is assumed to have the following columns:
#
#    time  xpos ypos zpos  xvel yvel zvel  xacc yacc zacc  Ekin/m Epot/m
#
# Nr:  1     2    3    4     5    6    7     8    9   10     11     12
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# plot #1: orbit

  r = 5
  set xrange [-r:r]
  set yrange [-r:r]
  set size ratio -1
  set xlabel 'x (internal units)'
  set ylabel 'y (internal units)'

  plot \
    'T3.out' u 2:3 w p lt 1 pt 7 ps .2 title 'orbit'

#-----------------------------------------------------------------------
# plot #2: energy vs. time

  set auto x
  set auto y
  set size noratio
  set xlabel 't (internal units)'
  set ylabel 'E/m (internal units)'

  plot \
    0 w l lt 2 lc rgb '#cccccc' notitle, \
    'T3.out' u 1:12        w l lt 3 title 'potential', \
    'T3.out' u 1:11        w l lt 1 title 'kinetic', \
    'T3.out' u 1:($11+$12) w l lt 7 title 'total'

#-----------------------------------------------------------------------

