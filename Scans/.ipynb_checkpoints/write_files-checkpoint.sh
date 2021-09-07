for angle in $(cd angle_xyzs; ls *.xyz; cd ../); do
    echo ${angle}

    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_template.com > AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.com
    sed "1d;2d" angle_xyzs/${angle} >> AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.com
    echo "" >> AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.com
    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_template.sub > AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.sub

    if [[ $angle == *"monomer"* ]]; then
        sed "s/NAME/${angle/.xyz/}/g" exciton_bchla_xtb_template.in > AngleScans/Exciton_Bchla_xTB/exciton_bchla_xtb_${angle/.xyz/}.in;
    fi

done
