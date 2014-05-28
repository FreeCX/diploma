set terminal pdf font 'Meslo LG L, 12' size 5.0in,4.0in
set format y '%.1f'
set key above
set output 'band_profile.pdf'
plot 'profileF1.dat' w l lt -1 lw 4 t '1', \
     'profileF2.dat' w l lt 9 lw 4 t '2',\
     'profileB.dat' w l lt 0 lw 5 t '3'