# Set terminal
set terminal epslatex color colortext dashed
#Strichstärken Inkscape
#0.15mm
#0.20mm
#0.30mm
#0.50mm
#0.70mm
#1.50mm
# width = 90mm
# Set size of plot area
set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.95
show margin
#set size 10/(5*2.54), 9/(3*2.54)
#  size: ^^ x in cm   ^^ y in cm (plot area only)
# set size square

#Set linestyles width = 2
set style line 1 lc rgb '#963821' pt 1 ps 1 lt 1 lw 2 # --- red  Pantone 174 C
set style line 2 lc rgb '#5e9c36' pt 2 ps 1 lt 2 lw 2 # --- green
set style line 3 lc rgb '#0065BD' pt 3 ps 1 lt 3 lw 2 # --- TUMblue
set style line 4 lc rgb '#FF8674' pt 4 ps 1 lt 4 lw 2 # --- orange Pantone 170 C
set style line 5 lc rgb '#A500CF' pt 5 ps 1 lt 5 lw 2 # --- magenta
set style line 6 lc rgb '#00998B' pt 6 ps 1 lt 6 lw 2 # --- cyan
set style line 7 lc rgb '#E2DA28' pt 7 ps 1 lt 7 lw 2 # --- yellow
set style line 8 lc rgb '#8A4B08' pt 8 ps 1 lt 8 lw 2 # --- brown
set style line 9 lc rgb '#088A29' pt 9 ps 1 lt 9 lw 2 # --- dark green
set style line 10 lc rgb '#FF0000' pt 10 ps 1 lt 10 lw 2 # --- red

# Set linestyles width = 2 and dashed
set style line 20 lc rgb '#963821' dt 2 pt 1 ps 1 lt 1 lw 2 # --- red  Pantone 174 C
set style line 21 lc rgb '#000000' dt 2 pt 1 ps 1 lt 1 lw 2 # --- black
set style line 22 lc rgb '#01DF01' dt 2 pt 2 ps 1 lt 2 lw 2 # --- light green
set style line 23 lc rgb '#0065BD' dt 2 pt 3 ps 1 lt 3 lw 2 # --- TUMblue
set style line 24 lc rgb '#FF8674' dt 2 pt 4 ps 1 lt 4 lw 2 # --- orange Pantone 170 C
set style line 25 lc rgb '#A500CF' dt 2 pt 5 ps 1 lt 5 lw 2 # --- magenta
set style line 26 lc rgb '#00998B' dt 2 pt 6 ps 1 lt 6 lw 2 # --- cyan
set style line 27 lc rgb '#E2DA28' dt 2 pt 7 ps 1 lt 7 lw 2 # --- yellow
set style line 28 lc rgb '#8A4B08' dt 2 pt 8 ps 1 lt 8 lw 2 # --- brown
set style line 29 lc rgb '#088A29' dt 2 pt 9 ps 1 lt 9 lw 2 # --- dark green
set style line 30 lc rgb '#FF0000' dt 2 pt 10 ps 1 lt 10 lw 2 # --- red
set style line 31 lc rgb '#5e9c36' dt 2 pt 1 ps 1 lt 1 lw 2 # # --- green

#Set linestyles width = 3
set style line 102 lc rgb '#963821' pt 1 ps 1 lt 1 lw 3 # --- red  Pantone 174 C
set style line 202 lc rgb '#5e9c36' pt 2 ps 1 lt 2 lw 3 # --- green
set style line 302 lc rgb '#0065BD' pt 3 ps 1 lt 3 lw 3 # --- TUMblue
set style line 402 lc rgb '#FF8674' pt 4 ps 1 lt 4 lw 3 # --- orange Pantone 170 C
set style line 502 lc rgb '#A500CF' pt 5 ps 1 lt 5 lw 3 # --- magenta
set style line 602 lc rgb '#00998B' pt 6 ps 1 lt 6 lw 3 # --- cyan
set style line 702 lc rgb '#E2DA28' pt 7 ps 1 lt 7 lw 3 # --- yellow
#Set linestyles width = 3 for matching pairs
set style line 1022 lc rgb '#963821' pt 2 ps 1 lt 1 lw 3 # --- red  Pantone 174 C
set style line 2022 lc rgb '#5e9c36' pt 6 ps 1 lt 2 lw 3 # --- green
set style line 3022 lc rgb '#0065BD' pt 2 ps 1 lt 3 lw 6 # --- TUMblue
set style line 4022 lc rgb '#FF8674' pt 6 ps 1 lt 4 lw 3 # --- orange Pantone 170 C
set style line 5022 lc rgb '#A500CF' pt 2 ps 1 lt 5 lw 3 # --- magenta
set style line 6022 lc rgb '#00998B' pt 6 ps 1 lt 6 lw 3 # --- cyan
set style line 7022 lc rgb '#E2DA28' pt 2 ps 1 lt 7 lw 3 # --- yellow
set style line 8022 lc rgb '#383838' pt 2 ps 1 lt 1 lw 3 # --- dark grey


#Tics and grid
set style line 11 lc rgb '#000000' lt 1 lw 2
set border 3 back ls 11
set tics nomirror
set style line 12 lc rgb '#808080' lt 0 lw 3
set grid back ls 12

#Key
#set key right top
unset key
unset title

