set terminal pngcairo size 900,650 enhanced font 'Verdana,10'
set output 'output/loglog_error_plot.png'
set title 'Log-Log Error Plot: Forward/Backward/Central Difference'
set xlabel 'log10(h)'
set ylabel 'log10(|error|)'
set grid
set key outside right
plot \
  'output/diff_data_exp.dat' using 1:2 with linespoints title 'exp forward', \
  'output/diff_data_exp.dat' using 1:3 with linespoints title 'exp backward', \
  'output/diff_data_exp.dat' using 1:4 with linespoints title 'exp central', \
  'output/diff_data_sin.dat' using 1:2 with linespoints title 'sin forward', \
  'output/diff_data_sin.dat' using 1:3 with linespoints title 'sin backward', \
  'output/diff_data_sin.dat' using 1:4 with linespoints title 'sin central', \
  'output/diff_data_poly.dat' using 1:2 with linespoints title 'poly forward', \
  'output/diff_data_poly.dat' using 1:3 with linespoints title 'poly backward', \
  'output/diff_data_poly.dat' using 1:4 with linespoints title 'poly central'
