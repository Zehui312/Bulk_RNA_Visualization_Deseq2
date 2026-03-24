#!/bin/bash
pwd=$(pwd)
metadata_file="${pwd}/meta_data.csv"

sample_name=$(csvcut -c sample_name ${metadata_file} | tail -n +2 | uniq)
output_dir=$(csvcut -c Output ${metadata_file} | tail -n +2 | uniq)

output_path=${output_dir}/$sample_name

if [ ! -d ${output_path} ]; then
    mkdir -p ${output_path}
fi
cd ${output_path}

#=================================================================
#+++++++++++++++++++++++Step 1 feature counting ++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/1_feature_count
cd ${output_path}/1_feature_count

bam_files=$(cat ${metadata_file} | tail -n +2 | cut -d "," -f 3| tr '\n' ' ')

minimum_fragment=$(csvcut -c minimum_fragment ${metadata_file} | tail -n +2 | uniq)
maximum_fragment=$(csvcut -c maximum_fragment ${metadata_file} | tail -n +2 | uniq)
count_features=$(csvcut -c count_features ${metadata_file} | tail -n +2 | uniq | tr -d '"')
gff_annotation=$(csvcut -c gff_annotation ${metadata_file} | tail -n +2 | uniq)
gff_file_path=$(csvcut -c GFF_file_path ${metadata_file} | tail -n +2 | uniq)

ln -s $bam_files .

input_bam=$(ls *.bam | tr '\n' ' ')
featureCounts -O -d ${minimum_fragment} -D ${maximum_fragment} -t ${count_features} -g ${gff_annotation} -a ${gff_file_path} -o ${sample_name}_gene.count -T 4 -L ${input_bam}


rm *.bam
grep -v "^#" ${sample_name}_gene.count > rename_${sample_name}_gene.count

csvcut -c "Bam_file_path","Sample_ID" ${metadata_file} | tail -n +2 | while IFS=, read bam_file sample_id; do
    filename=$(basename ${bam_file})
    sed -i "1s|${filename}|${sample_id}|g" rename_${sample_name}_gene.count
done

#=================================================================
#+++++++++++++++++++++++Step 2 Deseq2 +++++++++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/2_deseq2
cd ${output_path}/2_deseq2


deseq_rscript=${pwd}/script/Deseq2.R
cp ${deseq_rscript} .
cp ${output_path}/1_feature_count/rename_${sample_name}_gene.count .

Rscript Deseq2.R -c rename_${sample_name}_gene.count -m ${metadata_file} -o ${output_path}/2_deseq2
rm rename_${sample_name}_gene.count Deseq2.R Rplots.R