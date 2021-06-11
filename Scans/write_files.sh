for sep in $(cd sep_xyzs; ls *.xyz; cd ../); do
    echo ${sep}

    sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" ../templates/Bchla_template.in > DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.in
    sed "s/NAME/${sep/.xyz/}/g" ../templates/Bchla_template.sub > DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.sub
    sed -ie "s/DimerModel\/Bchla_xTB/DimerModel\/Scans\/DistanceScans\/Bchla_xTB/" DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.sub

    sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" ../templates/PBE0_template.in > DistanceScans/PBE0/PBE0_${sep/.xyz/}.in
    sed "s/NAME/${sep/.xyz/}/g" ../templates/PBE0_template.sub > DistanceScans/PBE0/PBE0_${sep/.xyz/}.sub
    sed -ie "s/DimerModel\/PBE0/DimerModel\/Scans\/DistanceScans\/PBE0/" DistanceScans/PBE0/PBE0_${sep/.xyz/}.sub
done

for angle in $(cd angle_xyzs; ls *.xyz; cd ../); do
    echo ${angle}

    sed "s/NAME/..\/..\/angle_xyzs\/${angle/.xyz/}/g" ../templates/Bchla_template.in > AngleScans/Bchla_xTB/Bchla_${angle/.xyz/}.in
    sed "s/NAME/${angle/.xyz/}/g" ../templates/Bchla_template.sub > AngleScans/Bchla_xTB/Bchla_${angle/.xyz/}.sub
    sed -ie "s/DimerModel\/Bchla_xTB/DimerModel\/Scans\/AngleScans\/Bchla_xTB/" AngleScans/Bchla_xTB/Bchla_${angle/.xyz/}.sub

    sed "s/NAME/..\/..\/angle_xyzs\/${angle/.xyz/}/g" ../templates/PBE0_template.in > AngleScans/PBE0/PBE0_${angle/.xyz/}.in
    sed "s/NAME/${angle/.xyz/}/g" ../templates/PBE0_template.sub > AngleScans/PBE0/PBE0_${angle/.xyz/}.sub
    sed -ie "s/DimerModel\/PBE0/DimerModel\/Scans\/AngleScans\/PBE0/" AngleScans/PBE0/PBE0_${angle/.xyz/}.sub
done