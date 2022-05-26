#for angle in $(cd angle_xyzs; ls *.xyz; cd ../); do

angles=(
sep_17_monomer_angle_0_axis_Qx.xyz
)
#sep_17_monomer_angle_0_axis_Qy.xyz
#sep_17_monomer_angle_0_axis_Qz.xyz
#sep_25_monomer_angle_0_axis_Qx.xyz
#sep_25_monomer_angle_0_axis_Qy.xyz
#sep_25_monomer_angle_0_axis_Qz.xyz
#sep_60_monomer_angle_0_axis_Qx.xyz
#sep_60_monomer_angle_0_axis_Qy.xyz
#sep_60_monomer_angle_0_axis_Qz.xyz
#sep_17_monomer_angle_181_axis_Qx.xyz
#sep_17_monomer_angle_181_axis_Qy.xyz
#sep_17_monomer_angle_181_axis_Qz.xyz
#sep_25_monomer_angle_181_axis_Qx.xyz
#sep_25_monomer_angle_181_axis_Qy.xyz
#sep_25_monomer_angle_181_axis_Qz.xyz
#sep_60_monomer_angle_181_axis_Qx.xyz
#sep_60_monomer_angle_181_axis_Qy.xyz
#sep_60_monomer_angle_181_axis_Qz.xyz
#)

for angle in ${angles[@]}; do
    echo ${angle}

    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_ground_template.com > AngleScans/CAMB3LYP/CAMB3LYP_ground_${angle/.xyz/}.com
    sed "1d;2d" angle_xyzs/${angle} >> AngleScans/CAMB3LYP/CAMB3LYP_ground_${angle/.xyz/}.com
    echo "" >> AngleScans/CAMB3LYP/CAMB3LYP_ground_${angle/.xyz/}.com
    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_ground_template.sub > AngleScans/CAMB3LYP/CAMB3LYP_ground_${angle/.xyz/}.sub

    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_transition_template.com > AngleScans/CAMB3LYP/CAMB3LYP_transition_${angle/.xyz/}.com
    sed "1d;2d" angle_xyzs/${angle} >> AngleScans/CAMB3LYP/CAMB3LYP_transition_${angle/.xyz/}.com
    echo "" >> AngleScans/CAMB3LYP/CAMB3LYP_transition_${angle/.xyz/}.com
    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_transition_template.sub > AngleScans/CAMB3LYP/CAMB3LYP_transition_${angle/.xyz/}.sub

    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_excited_template.com > AngleScans/CAMB3LYP/CAMB3LYP_excited_${angle/.xyz/}.com
    sed "1d;2d" angle_xyzs/${angle} >> AngleScans/CAMB3LYP/CAMB3LYP_excited_${angle/.xyz/}.com
    echo "" >> AngleScans/CAMB3LYP/CAMB3LYP_excited_${angle/.xyz/}.com
    sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_excited_template.sub > AngleScans/CAMB3LYP/CAMB3LYP_excited_${angle/.xyz/}.sub

    #if [[ $angle == *"monomer"* ]]; then
    #    sed "s/NAME/${angle/.xyz/}/g" exciton_bchla_xtb_template.in > AngleScans/Exciton_Bchla_xTB/exciton_bchla_xtb_${angle/.xyz/}.in;
    #fi

done
