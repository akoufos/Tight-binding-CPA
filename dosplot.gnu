#    G N U P L O T
#    Linux version 3.5 (pre 3.6)
#    patchlevel beta 347
#    last modified Mon Jun 22 13:22:33 BST 1998
#
#    Copyright(C) 1986 - 1993, 1998
#    Thomas Williams, Colin Kelley and many others
#
#    Send comments and requests for help to info-gnuplot@dartmouth.edu
#    Send bugs, suggestions and mods to bug-gnuplot@dartmouth.edu
# Name: dosplot.gnu
# Objective: plot density of states of monoatomic materials
# Alex Koufos, 05/16/2013
fermi = 0.30500 
set xlabel "Energy(Ry)"
set ylabel "States/Ry/Atom"

set style line 1 lt 1 lw 1 lc 1
set style line 2 lt 2 lw 1 lc 2
set style line 3 lt 3 lw 1 lc 3
set style line 4 lt 4 lw 1 lc 4
set style line 5 lt 5 lw 1 lc 5
set style line 6 lt 2 lw 1 lc 6
set style line 7 lt 3 lw 1 lc 7
set style line 8 lt 4 lw 1 lc 8
set style line 9 lt 5 lw 1 lc 9

set key left top Right noreverse box linetype -2 linewidth 1.000 samplen 4 spacing 1 width 0
set nolabel
set noarrow
set nologscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default

#Plot total DOS and all angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):10 t "Total DOS" w l ls 1,\
     "dosdat.cpa.plot" u ($1-fermi):2 t "DOS Fe-s" w l ls 2,\
     "dosdat.cpa.plot" u ($1-fermi):3 t "DOS Fe-p" w l ls 3,\
     "dosdat.cpa.plot" u ($1-fermi):($4+$5+$6+$7) t "DOS Fe-d" w l ls 4, \
     "dosdat.cpa.plot" u ($1-fermi):8 t "DOS Se-s" w l ls 5, \
     "dosdat.cpa.plot" u ($1-fermi):9 t "DOS Se-p" w l ls 6

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatall.eps"
replot

#Plot total DOS 
plot "dosdat.cpa.plot" u ($1-fermi):10 t "Total DOS" w l ls 1, \
     "dosdat.cpa.plot" u ($1-fermi):($2+$8) t "DOS total-s" w l ls 2,\
     "dosdat.cpa.plot" u ($1-fermi):($3+$9) t "DOS total-p" w l ls 3,\
     "dosdat.cpa.plot" u ($1-fermi):($4+$5+$6+$7) t "DOS total-d" w l ls 4

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdattotal.eps"
replot

#Plot total angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):($2+$8) t "DOS total-s" w l ls 2,\
     "dosdat.cpa.plot" u ($1-fermi):($3+$9) t "DOS total-p" w l ls 3,\
     "dosdat.cpa.plot" u ($1-fermi):($4+$5+$6+$7) t "DOS total-d" w l ls 4

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatang.eps"
replot

#Plot total S+P angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):($2+$8) t "DOS total-s" w l ls 2,\
     "dosdat.cpa.plot" u ($1-fermi):($3+$9) t "DOS total-p" w l ls 3

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatsp.eps"
replot

#Plot total S angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):($2+$8) t "DOS total-s" w l ls 2

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdats.eps"
replot

#Plot total P angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):($3+$9) t "DOS total-p" w l ls 3

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatp.eps"
replot

#Plot total D angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):($4+$5+$6+$7) t "DOS total-d" w l ls 4

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdatd.eps"
replot

#Plot total DOS and atom kind 1's angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):10 t "Total DOS" w l ls 1,\
     "dosdat.cpa.plot" u ($1-fermi):2 t "DOS Fe-s" w l ls 2,\
     "dosdat.cpa.plot" u ($1-fermi):3 t "DOS Fe-p" w l ls 3,\
     "dosdat.cpa.plot" u ($1-fermi):($4+$5+$6+$7) t "DOS Fe-d" w l ls 4

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdat1all.eps"
replot

#Plot total DOS and atom kind 2's angular momentum decompositions
plot "dosdat.cpa.plot" u ($1-fermi):10 t "Total DOS" w l ls 1,\
     "dosdat.cpa.plot" u ($1-fermi):2 t "DOS Se-s" w l ls 6,\
     "dosdat.cpa.plot" u ($1-fermi):3 t "DOS Se-p" w l ls 7

set label 1 "{/Symbol e}_F" at 0, 9*GPVAL_DATA_Y_MAX/10, 0 right norotate
set arrow 1 from 0, graph 0, 0 to 0, graph 1, 0  nohead linetype 3 linewidth 1.000
set term postscript eps enhanced color
set output "dosdat2all.eps"
replot

set output      
set terminal x11        
#    EOF
