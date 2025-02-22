set term gif animate delay 5 background "black"
set output "test.gif"

L = 3000
set xrange [-L:L]
set yrange [-L:L]
unset colorbox
unset xlabel
unset ylabel
unset xtics
unset ytics
set border linecolor "white"
set size ratio -1

n = 10000
T = 10000
  
do for [i=0:T*n:100*n] {
    #filename = sprintf("quadranti/quadrants_%d.dat", i/n)
    plot "dati.dat" u 1:2 every ::i::i+(n-1) pt 7 ps 0.1 lc "white" notitle, \
         #filename using 1:2 with lines lc "green" notitle
}