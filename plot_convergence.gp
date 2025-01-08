# Set terminal and output
set terminal pngcairo enhanced font 'Arial,12' size 800,600
set output 'convergence_comparison.png'

# Set title and labels
set title "Convergence History Comparison"
set xlabel "Iteration"
set ylabel "Residual"

# Set grid and style
set grid

# Set ranges
set xrange [0:140]
set yrange [0:1]

# Set line styles
set style line 1 lc rgb '#0000FF' lt 1 lw 2 pt 7 ps 0.5  # Blue for Richardson
set style line 2 lc rgb '#FF0000' lt 1 lw 2 pt 9 ps 0.5  # Red for Jacobi
set style line 3 lc rgb '#00AA00' lt 1 lw 2 pt 5 ps 0.5  # Green for Gauss-Seidel

# Plot all three methods
plot 'RESVEC_alpha.dat' using 0:1 with linespoints ls 1 title 'Richardson', \
     'RESVEC_jacobi.dat' using 0:1 with linespoints ls 2 title 'Jacobi', \
     'RESVEC_GS.dat' using 0:1 with linespoints ls 3 title 'Gauss-Seidel'