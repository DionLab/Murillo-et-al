#BWA

for fastq_r1 in "${input_dir}"*_R1_001.fastq.gz; do
    fastq_r2="${fastq_r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    base_name=$(basename "${fastq_r1}" _R1_001.fastq.gz)

    output_sam="${output_dir}/${base_name}.sam"
    output_bam="${output_dir}/${base_name}.sorted.bam"

    echo "[Mapping running for] $base_name" >> "${output_dir}/map_pairs.log"

    # BWA MEM Alignment
    bwa mem -M -t 36 "$REFERENCE_GENOME" "$fastq_r1" "$fastq_r2" > "$output_sam"
    echo "BWA MEM completed for $base_name" >> "${output_dir}/map_pairs.log"

    # Use Picard to convert SAM to BAM, sort and clean up
    java -jar $PICARD SortSam \
        INPUT="$output_sam" \
        OUTPUT="$output_bam" \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_RECORDS_IN_RAM=250000 \
        TMP_DIR="$tmp_dir"

    # Remove the temporary SAM file
    rm "$output_sam"

    echo "Processing completed for $base_name" >> "${output_dir}/map_pairs.log"
done


#Adding Read group

# Iterate over each .bam file in the source directory
for file in "$SOURCE_DIR"/*.bam; do
    # Extract the base name without the .sort_rg.bam extension
    base_name=$(basename "$file" .sort_rg.bam)

    # Check if sort_rg.bam exists for the current sample
    if [ -f "$DESTINATION_DIR/${base_name}.sort_rg.bam" ]; then
        echo "sort_rg.bam exists for ${base_name}, skipping processing."
    else
        echo "Processing: $base_name"

        # Run the PICARD command to add or replace read groups
        java -jar "$PICARD" AddOrReplaceReadGroups \
        I="${file}" \
        O="${DESTINATION_DIR}/${base_name}.sort_rg.bam" \
        RGID=4 \
        RGLB=lib1 \
        RGPL=ILLUMINA \
        RGPU=unit1 \
        RGSM=20

        echo "Processed: ${base_name}.sort_rg.bam"
    fi
done

echo "Processing completed."

#Indexing bam file

for bam_file in "$bam_dir"/*.sort_rg.bam; do
    bai_file="${bam_file}.bai"
    if [ -s "$bam_file" ]; then
        if [ -e "$bai_file" ]; then
            echo "Skipping $bam_file as $bai_file already exists."
        else
            echo "Indexing $bam_file..."
            samtools index "$bam_file"
        fi
    else
        echo "Skipping empty or non-existent file: $bam_file"
    fi
done

#Mosdepth for coverage

for i in /scratch/scw1617/Ruban_backup_2024/S4_WGS_dec23/read_assigned_files/*.bam; do
    base=$(basename "$i" .bam)
    output_path="${output_directory}/${base}-interval"

    # Run Mosdepth
    mosdepth -n --by $bed_file --fast-mode $output_path $i

    # Unzip and save the regions file
    gzip -dc "${output_path}.regions.bed.gz" > "${output_path}.regions.bed"
done


#Variant Calling

# Iterate over each BAM file in the specified directory
for BAM_FILE in /scratch/scw1617/Ruban_backup_2024/S4_WGS_dec23/bam_files/*.bam; do
    # Extract the base name of the BAM file (without .bam extension)
    BAM_BASENAME=$(basename "${BAM_FILE}" .bam)

    # Create a dedicated output directory for each BAM file
    OUTPUT_DIR="${OUTPUT_DIR_BASE}/${BAM_BASENAME}"

    # If the output directory already exists, skip processing this BAM file
    if [ -d "${OUTPUT_DIR}" ]; then
        echo "Output folder ${OUTPUT_DIR} exists. Skipping ${BAM_BASENAME}."
        continue
    fi

    # Make output directory
    mkdir -p "${OUTPUT_DIR}"

    # Define output file paths
    RAW_BCF="${OUTPUT_DIR}/raw_calls.bcf"
    FILTERED_VCF="${OUTPUT_DIR}/filtered_calls.vcf.gz"

    # Step 1: Generate genotype likelihoods and call variants in specified regions
    bcftools mpileup -Ou -f "${REFERENCE}" "${BAM_FILE}" | \
        bcftools call -mv -Ob -o "${RAW_BCF}"

    # Step 2: Filter variants based on quality and other criteria
    bcftools filter -e 'QUAL<30 || DP<10 || (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>=0.2 || MQ<40' "${RAW_BCF}" -Ou | \
        bcftools view -e 'GT="0/1"' -Oz -o "${FILTERED_VCF}"

    # Step 3: Index the filtered VCF
    bcftools index "${FILTERED_VCF}"

    # Step 4: Generate statistics for the filtered VCF
    bcftools stats "${FILTERED_VCF}" > "${OUTPUT_DIR}/filtered_stats.txt"

    echo "Variant calling and filtering completed for ${BAM_FILE}. Output files saved in ${OUTPUT_DIR}"
done