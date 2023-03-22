#!/usr/bin/bash
# Grab croco*.out files and extract kinetic, potential and total energy time series on execution time.


yrstart=2000
yrend=2001
mstart=1
mend=12
sim_name='crocod1'
outdirectory="/ceaza/lucas/CROCO/DESALADORAS_RUND1/OUTPUT"
for yr in $(seq $yrstart $yrend); do
	for m in $(seq ${mstart} ${mend}); do
		timestamp=Y${yr}M${m}
		f=$outdirectory/${sim_name}_${timestamp}.out
		echo $f
		START=$(cat $f | grep -n 'STEP' | awk -F ':' '{print $1}')
		TLINES=$(cat $f | wc -l)
		START=$(( $TLINES - $START + 1 ))
		tail -n $START $f > .tmp
		END=$(cat .tmp | grep -n 'DEF_RST' | awk -F ':' '{print $1}')
		cat .tmp | head -n $(( $END - 1 )) > .${sim_name}_status.${timestamp}
		sed -i '/BRY/d' .${sim_name}_status.${timestamp}
		sed -i '/BULK/d' .${sim_name}_status.${timestamp}
		sed -i '/AVG/d' .${sim_name}_status.${timestamp}
		sed -i '/HIS/d' .${sim_name}_status.${timestamp}
		sed -i '/WRT/d' .${sim_name}_status.${timestamp}
		sed -i '/dapl/d' .${sim_name}_status.${timestamp}
		sed -i 's/  //g' .${sim_name}_status.${timestamp}
	done
done


head -n .${sim_name}_status_Y${yrstart}M${mstart} > ${sim_name}_status.txt
for f in .${sim_name}*; do
	tail -n +1 $f >> ${sim_name}_status.out
done

rm .${sim_name}*		


