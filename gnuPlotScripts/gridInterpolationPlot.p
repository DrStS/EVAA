cd "c:/software/repos/EVAA/build/src/"
set datafile separator ","
set key width 6
set key height 1
set key box
set key font ",18"
set ylabel "Stiffness in [Nm]" font ",18"
set xlabel "Spring length in [m]" font ",18"
plot "LookupTablePlusInterpolation.txt" us 1:2 lt 7 lc 7 w lines title "interpolation", "LookupTablePlusInterpolation.txt" using 3:4 lt 1 lc 1 title "grid"