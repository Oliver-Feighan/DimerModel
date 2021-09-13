for frame in $(seq 1 100 49951); do
    for i in $(seq 1 27); do
        for j in $(seq ${i} 27); do
        
            if [ ${i} == ${j} ] ; then
                continue
            fi;
            
            full_name=trunc_bchla_${i}_bchla_${j}_frame_${frame}
            
            echo $full_name
            
            monomerA=trunc_bchla_${i}_frame_${frame}.xyz
            monomerB=trunc_bchla_${j}_frame_${frame}.xyz
            
            sed -e "s/FULL_NAME/${full_name}/g" -e "s/NAMEA/${monomerA}/g" -e "s/NAMEB/${monomerB}/g" templates/exciton_bchla_xtb_template.in > Exciton_Bchla_xTB/exciton_bchla_xtb_${full_name}.in
            
        done;
    done;
done;
    
        