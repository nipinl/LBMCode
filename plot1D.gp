set terminal postscript eps color enhanced
set terminal jpeg enhanced size 800,600
set output 'T.jpg';
set xlabel "Length"font "Times-Roman,20";
set ylabel "Temperature"font "Times-Roman,20";
set title "D1Q3" font "Times-Roman,24";
plot "./T"  using 1:2 t 'Temperature'  w lp pt 7 lt 7 lw 4
