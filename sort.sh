#!/bin/bash
#sort.sh: bash script to sort the LE4PD modes

for j in 1 2 3 
do

	for i in 0.1 0.2 0.4 0.5 0.6 0.8 1.0 2.0
	do

		cd 1UBQRun${j}/1UBQ${i}exp/
		echo "Sorting run ${j} and relative viscosity ${i}"
		cp -v ../../sort.py ./
		python sort.py
		nlines=$( wc -l tau_sorted.dat  | awk '{ print $1 }' )
		n=1
		while read line
		do
			echo "$line" >line
			mode=$( awk '{ print $1 }' line)
			echo "$mode $n" >> holder.dat
			n=$(( $n + 1 ))
		done < "tau_sorted.dat"	
		f95 ../../mode_sort.f95
		./a.out
		#Re-write free-energy surfaces:
		while read line
		do
			echo "$line" >line
			oldmode=$( awk '{ print $1 }' line)
			newmode=$( awk '{ print $2 }' line)
			cp fempa_${oldmode}.dat new_fempa_${newmode}.dat
			cp fe_${oldmode}.dat new_fe_${newmode}.dat
		done < "old_mode_new_mode.dat"	
		#Re-write barriers:
		f95 ../../barrier_sort.f95
		./a.out		
		cd ../..
	done
done

