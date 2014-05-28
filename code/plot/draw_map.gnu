set terminal png font 'Meslo LG L' 12 size 800,800
set size square 1.0,1.0
set palette rgbformulae 33,13,10
set xrange[300:900]
set yrange[300:900]
set format z '%.1f'
set format cb '%.1f'
# set contour base
# set cntrparam level auto 10
# set noclabel
set notitle
set pm3d map

set output 'map_B.png'
splot 'dataB.dat' every 1:1 matrix w pm3d palette title ''

set output 'map_F1.png'
splot 'dataF1.dat' every 1:1 matrix w pm3d palette title ''

set output 'map_F2.png'
splot 'dataF2.dat' every 1:1 matrix w pm3d palette title ''