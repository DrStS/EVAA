load "../common_template.gp"

set output "insulationZeller2Dof.eps"
set xlabel 'Stiffness Ratio'
set ylabel 'Difference level D [dB]'
set xrange [0:20]
set xrange [0:30]

plot 'data/2dofZeller.csv' ls 1 with lines title '$ 20 \cdot \log \left| \nicefrac{v_o}{v_u} \right| $', \
		 'data/2dofZellerAnalytic.csv' ls 22 with lines title '$ 20 \cdot \log \left| \nicefrac{\left(k_1 + k_2\right)}{k_1}\right| $'