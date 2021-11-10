# Gnuplot script file
set terminal postscript eps color enhanced "Arial" 20 size 7in,7in
set output "radiation_force2.eps"

#set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set pm3d map
set multiplot layout 2,1

set grid lc rgb "white" 
set xrange[10e3 : 200e3]
set yrange[1e-3 : 20e-3]

set title "z-component of radiation force"
set xlabel "frequency [Hz]"
set ylabel "sphere radius [m]"
set cblabel "radiation force [N]"
splot "radiation_force2.txt"

set title "plot only negative value (attractive force toward the wave source)"
set cbrange[ : 0]
set palette negative
splot "radiation_force2.txt"
