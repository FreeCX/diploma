set terminal pdf font 'Meslo LG L, 12' size 5.0in,4.0in
set format y '%.1f'
set key above
set output 'band_profile.pdf'
plot 'profileF1.dat' w l lw 3 t 'зона 1', \
     'profileF2.dat' w l lw 3 t 'зона 2',\
     'profileB.dat' w l lw 3 t 'поле B'