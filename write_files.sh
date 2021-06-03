for frame in $(seq 1 50 151); do
	echo ${frame}
	
	for i in $(seq 1 27); do
		
		monomer=trunc_bchla_${i}_frame_${frame}
		echo $monomer

		sed "s/NAME/${monomer}/g" templates/BCHLA_template.in > Bchla_xTB/Bchla_${monomer}.in
		sed "s/NAME/${monomer}/g" templates/BCHLA_template.sub > Bchla_xTB/Bchla_${monomer}.sub


		for j in $(seq ${i} 27); do
			
			if [ ${i} == ${j} ] ; then
				continue
			fi;

			dimer=trunc_bchla_${i}_bchla_${j}_frame_${frame}
			echo $dimer

			sed "s/NAME/${dimer}/g" templates/PBE0_template.in > PBE0/PBE0_${dimer}.in	
			sed "s/NAME/${dimer}/g" templates/PBE0_template.sub > PBE0/PBE0_${dimer}.sub


			sed "s/NAME/${dimer}/g" templates/Bchla_template.in > Bchla_xTB/Bchla_${dimer}.in	
			sed "s/NAME/${dimer}/g" templates/Bchla_template.sub > Bchla_xTB/Bchla_${dimer}.sub
		done
	done
done		
