unset output

# set colors classic
set datafile separator ' '

scale = 3.5
rows = 999
set xrange [-scale:scale]
set yrange [-scale:scale]
set size ratio 1
set xlabel "x (AU)"
set ylabel "y (AU)"
set grid

plot "dt1e-2-position-first.txt" u 2:3 every ::0::rows title "Earth Orbit, dt=1e-2 (year)" w p pt 7 ps .5 lc rgb 'blue',\
     "center.txt" u 1:2 title "the Sun" w p pt 7 ps 3 lc rgb 'red', \
    #  "dt1e-2-position-first.txt" u 2:3 every ::0::rows title "Earth Orbit, dt=1e-2 (year)" w l dt 2 lw 2 lc rgb 'grey'