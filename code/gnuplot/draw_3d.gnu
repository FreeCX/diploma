set terminal png font 'Meslo LG L' 12 size 800,800
set size square 1.0,1.0
set palette rgbformulae 33,13,10
set xrange[0:1200]
set yrange[0:1200]
set format z '%.1f'
set format cb '%.1f'
set ytic offset character +2,+0
set xtic offset character -2,-1
set ticslevel 0
set noclabel
set pm3d

set output '3d_B.png'
set title 'поле B'
splot 'dataB.dat' every 1:1 matrix w pm3d palette title ''

set contour base
set cntrparam level auto 10

set output '3d_F1.png'
set title 'зона 1'
splot 'dataF1.dat' every 1:1 matrix w pm3d palette title ''

set output '3d_F2.png'
set title 'зона 2'
splot 'dataF2.dat' every 1:1 matrix w pm3d palette title ''