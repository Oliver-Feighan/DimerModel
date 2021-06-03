#run in this clean_pdbs dir

for i in ../raw_pdbs/*; do
	echo ${i}
	clean_file=${i/..\/raw_pdbs\//clean_}
	clean_file=${clean_file/.pdb./_frame_}
	sed "s/ '/  /g" ${i} > ${clean_file}.pdb 
done
