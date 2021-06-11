for sep in $(cd sep_xyzs; ls *.xyz; cd ../); do
    echo ${sep}

    sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" ../templates/Bchla_template.in > DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.in
    sed "s/NAME/${sep/.xyz/}/g" ../templates/Bchla_template.sub > DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.sub
    sed -ie "s/DimerModel\/Bchla_xTB/DimerModel\/Scans\/DistanceScans\/Bchla_xTB/" DistanceScans/Bchla_xTB/Bchla_${sep/.xyz/}.sub

    sed "s/NAME/..\/..\/sep_xyzs\/${sep/.xyz/}/g" ../templates/PBE0_template.in > DistanceScans/PBE0/PBE0_${sep/.xyz/}.in
    sed "s/NAME/${sep/.xyz/}/g" ../templates/PBE0_template.sub > DistanceScans/PBE0/PBE0_${sep/.xyz/}.sub
    sed -ie "s/DimerModel\/PBE0/DimerModel\/Scans\/DistanceScans\/PBE0/" DistanceScans/PBE0/PBE0_${sep/.xyz/}.sub
done
