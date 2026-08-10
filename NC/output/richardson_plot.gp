set terminal pngcairo size 1000,700 enhanced font 'Verdana,10'
set output 'output/richardson_loglog_plot.png'
set title 'Log-Log Error Plot: Central Difference D(h) vs Richardson R(h)'
set xlabel 'log10(h)'
set ylabel 'log10(|error|)'
set grid
set key outside right
plot \
  'output/richardson_data_sin.dat' using 1:2 with linespoints title 'sin D(h) [O(h^2)]', \
  'output/richardson_data_sin.dat' using 1:3 with linespoints title 'sin R(h) [O(h^4)]'
