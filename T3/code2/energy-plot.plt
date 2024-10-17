unset output

set datafile separator ' '

set xrange [0:1]
# set yrange [-0.0002:0]

set title "Energy Of Earth"
set xlabel "time (year)"
set ylabel "energy (M_sun * AU^2 / year^2)"
set grid

# plot "dt1e-1-with-E.txt" u 1:5 every ::0::9 title "dt=1e-1 (year)" w l lw 2 lc rgb 'red',\
#      "dt1e-2-with-E.txt" u 1:5 every ::0::99 title "dt=1e-2 (year)" w l lw 2 lc rgb 'blue',\
#      "dt1e-3-with-E.txt" u 1:5 every ::0::999 title "dt=1e-3 (year)" w l lw 2 lc rgb 'green'
plot "dt1e-1-position-first.txt" u 1:5 every ::0::9 title "dt=1e-1 (year)" w l lw 2 lc rgb 'red',\
     "dt1e-2-position-first.txt" u 1:5 every ::0::99 title "dt=1e-2 (year)" w l lw 2 lc rgb 'blue',\
     "dt1e-3-position-first.txt" u 1:5 every ::0::999 title "dt=1e-3 (year)" w l lw 2 lc rgb 'green'
