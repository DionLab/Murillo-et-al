#!/bin/bash

# Define input and output directories
variant_dir="/scratch/c.mpmrrp/WGS_projects/S4_WGS_dec23/bin/output/WG_variant_calling/"
output_csv="detailed_variant_summary.csv"

# Initialize the CSV file with headers
echo "Sample,#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,GENOTYPE,VAF" > "$output_csv"

# Loop through all sample directories
for sample_dir in "$variant_dir"*/; do
    sample_name=$(basename "$sample_dir")
    vcf_file="$sample_dir/filtered_calls.vcf.gz"

    # Check if the VCF file exists
    if [[ -f "$vcf_file" ]]; then
        # Extract detailed variant information and calculate VAF
        zgrep -v "^#" "$vcf_file" | awk -v sample="$sample_name" -F'\t' '
        BEGIN { OFS="," }
        {
            chrom = $1; pos = $2; id = $3; ref = $4; alt = $5; qual = $6;
            filter = $7; info = $8; format = $9; genotype = $10;

            # Parse DP4 to calculate VAF
            if (match(info, /DP4=([0-9]+),([0-9]+),([0-9]+),([0-9]+)/, dp4)) {
                ref_count = dp4[1] + dp4[2];
                alt_count = dp4[3] + dp4[4];
                total_count = ref_count + alt_count;
                if (total_count > 0) {
                    vaf = alt_count / total_count;
                } else {
                    vaf = "NA";
                }
            } else {
                vaf = "NA";
            }

            # Output details
            print sample, chrom, pos, id, ref, alt, qual, filter, info, format, genotype, vaf;
        }' >> "$output_csv"
    else
        # If the VCF file doesn't exist, add an entry with "No VCF"
        echo "$sample_name,,,,,,,,,,,No VCF" >> "$output_csv"
    fi
done

echo "Detailed variant summary written to $output_csv"

