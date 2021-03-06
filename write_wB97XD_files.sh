for frame in $(seq 1 250 501); do
    for i in $(seq 1 27); do

	monomer=trunc_bchla_${i}_frame_${frame}
	com_file=wB97XD/wB97XD_${monomer}.com
	sub_file=wB97XD/wB97XD_${monomer}.sub

	sed "s/NAME/${monomer}/g" templates/wB97XD_template.com > $com_file
	sed "1d;2d" monomer_xyzs/${monomer}.xyz >> $com_file
	echo >> $com_file
	
	sed "s/NAME/${monomer}/g" templates/wB97XD_template.sub > $sub_file

        for j in $(seq ${i} 27); do
            
            if [ ${i} == ${j} ] ; then
                continue
            fi;

            dimer=trunc_bchla_${i}_bchla_${j}_frame_${frame}
            echo $dimer

            com_file=wB97XD/wB97XD_${dimer}.com
            sub_file=wB97XD/wB97XD_${dimer}.sub

            sed "s/NAME/${dimer}/g" templates/wB97XD_template.com > $com_file
            sed "1d;2d" dimer_xyzs/${dimer}.xyz >> $com_file
            echo >> $com_file

            sed "s/NAME/${dimer}/g" templates/wB97XD_template.sub > $sub_file

        done
    done
done
