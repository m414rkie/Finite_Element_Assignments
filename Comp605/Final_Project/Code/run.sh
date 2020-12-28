#! /bin/bash

for i in 2 4 6 8 10 12
do
	eval "mpirun -np $i ./p2sparr.x < infile"
	eval "mv outs/coralperc.dat outs/coralperc$i.dat"
	eval "mv outs/fishtime.dat outs/fishtime$i.dat"
done

