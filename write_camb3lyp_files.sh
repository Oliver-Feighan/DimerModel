for frame in $(seq 1 250 501); do
    for i in $(seq 1 27); do

	monomer=trunc_bchla_${i}_frame_${frame}
	ground_com_file=CAMB3LYP/CAMB3LYP_ground_${monomer}.com
	ground_sub_file=CAMB3LYP/CAMB3LYP_ground_${monomer}.sub
	
	transition_com_file=CAMB3LYP/CAMB3LYP_transition_${monomer}.com
	transition_sub_file=CAMB3LYP/CAMB3LYP_transition_${monomer}.sub
	
	excited_com_file=CAMB3LYP/CAMB3LYP_excited_${monomer}.com
	excited_sub_file=CAMB3LYP/CAMB3LYP_excited_${monomer}.sub
	
	sed "s/NAME/ground_${monomer}/g" templates/CAMB3LYP_ground_template.com > $ground_com_file
	sed "1d;2d" monomer_xyzs/${monomer}.xyz >> $ground_com_file
	echo >> $ground_com_file
	
	sed "s/NAME/ground_${monomer}/g" templates/CAMB3LYP_ground_template.sub > $ground_sub_file


	sed "s/NAME/transition_${monomer}/g" templates/CAMB3LYP_template.com > $transition_com_file
	sed "1d;2d" monomer_xyzs/${monomer}.xyz >> $transition_com_file
	echo >> $transition_com_file
	
	sed "s/NAME/transition_${monomer}/g" templates/CAMB3LYP_transition_template.sub > $transition_sub_file


	sed "s/NAME/excited_${monomer}/g" templates/CAMB3LYP_excited_template.com > $excited_com_file
	sed "1d;2d" monomer_xyzs/${monomer}.xyz >> $excited_com_file
	echo >> $excited_com_file
	
	sed "s/NAME/excited_${monomer}/g" templates/CAMB3LYP_excited_template.sub > $excited_sub_file


        for j in $(seq ${i} 27); do
            
            if [ ${i} == ${j} ] ; then
                continue
            fi;

            #dimer=trunc_bchla_${i}_bchla_${j}_frame_${frame}
            #echo $dimer

            #com_file=CAMB3LYP/CAMB3LYP_${dimer}.com
            #sub_file=CAMB3LYP/CAMB3LYP_${dimer}.sub

            #sed "s/NAME/${dimer}/g" templates/CAMB3LYP_template.com > $com_file
            #sed "1d;2d" dimer_xyzs/${dimer}.xyz >> $com_file
            #echo >> $com_file

            #sed "s/NAME/${dimer}/g" templates/CAMB3LYP_template.sub > $sub_file

        done
    done
done
