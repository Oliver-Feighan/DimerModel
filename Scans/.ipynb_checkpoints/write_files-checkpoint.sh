for sep in $(cd sep_xyzs; ls *.xyz; cd ../); do
    echo ${sep}

    #sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" ../templates/Bchla_template.in > DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.in
    #sed "s/NAME/${sep/.xyz/}/g" ../templates/Bchla_template.sub > DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.sub
    #sed -ie "s/DimerModel\/Bchla_xTB/DimerModel\/Scans\/DistanceScans\/Bchla_xTB/" DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.sub

    #sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" CAMB3LYP_template.in > DistanceScans/CAMB3LYP/CAMB3LYP_${sep/.xyz/}.in
    #sed "s/NAME/${sep/.xyz/}/g" CAMB3LYP_template.sub > DistanceScans/CAMB3LYP/CAMB3LYP_${sep/.xyz/}.sub
    #sed -i "s/DimerModel\/CAMB3LYP/DimerModel\/Scans\/DistanceScans\/CAMB3LYP/" DistanceScans/CAMB3LYP/CAMB3LYP_${sep/.xyz/}.sub

    #sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" BLYP_template.in > DistanceScans/BLYP/BLYP_${sep/.xyz/}.in
    #sed "s/NAME/${sep/.xyz/}/g" BLYP_template.sub > DistanceScans/BLYP/BLYP_${sep/.xyz/}.sub
    #sed -i "s/DimerModel\/BLYP/DimerModel\/Scans\/DistanceScans\/BLYP/" DistanceScans/BLYP/BLYP_${sep/.xyz/}.sub

done

for angle in $(cd angle_xyzs; ls *.xyz; cd ../); do
    echo ${angle}

    sed "s/NAME/..\/..\/angle_xyzs\/${angle/.xyz/}/g" Bchla_template.in > AngleScans/Bchla_xTB/Bchla_${angle/.xyz/}.in
    sed "s/NAME/${angle/.xyz/}/g" Bchla_template.sub > AngleScans/Bchla_xTB/Bchla_${angle/.xyz/}.sub
    sed -ie "s/DimerModel\/Bchla_xTB/DimerModel\/Scans\/AngleScans\/Bchla_xTB/" AngleScans/Bchla_xTB/Bchla_${angle/.xyz/}.sub

    #sed "s/NAME/..\/..\/angle_xyzs\/${angle/.xyz/}/g" CAMB3LYP_template.in > AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.in
    #sed "s/NAME/${angle/.xyz/}/g" CAMB3LYP_template.sub > AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.sub
    #sed -i "s/DimerModel\/CAMB3LYP/DimerModel\/Scans\/AngleScans\/CAMB3LYP/" AngleScans/CAMB3LYP/CAMB3LYP_${angle/.xyz/}.sub

    #sed "s/NAME/..\/..\/angle_xyzs\/${angle/.xyz/}/g" BLYP_template.in > AngleScans/BLYP/BLYP_${angle/.xyz/}.in
    #sed "s/NAME/${angle/.xyz/}/g" BLYP_template.sub > AngleScans/BLYP/BLYP_${angle/.xyz/}.sub
    #sed -i "s/DimerModel\/BLYP/DimerModel\/Scans\/AngleScans\/BLYP/" AngleScans/BLYP/BLYP_${angle/.xyz/}.sub

    sed "s/NAME/..\/..\/angle_xyzs\/${angle/.xyz/}/g" HF_template.in > AngleScans/HF/HF_${angle/.xyz/}.in
    sed "s/NAME/${angle/.xyz/}/g" HF_template.sub > AngleScans/HF/HF_${angle/.xyz/}.sub
    sed -i -e "s/DimerModel\/HF/DimerModel\/Scans\/AngleScans\/HF/" AngleScans/HF/HF_${angle/.xyz/}.sub

done
