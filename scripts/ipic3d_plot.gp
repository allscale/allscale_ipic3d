set title "Simple demo of scatter data conversion to grid data"
set terminal windows size 1280,1024 title "AllScale ipic3d initial particle distribution"
unset hidden3d
set ticslevel 0.5
set view 60,30
set autoscale
set parametric
set style data points
set xlabel "data style point - no dgrid"
set key box
splot filename pt 7 ps 3
pause -1