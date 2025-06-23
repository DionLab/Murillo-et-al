import pandas as pd
import re

def reorder_chromosome_matrix(input_file, output_file):
    """
    Reorder chromosome coverage matrix columns by sample number and lane
    """
    
    # Read the original matrix
    print("Reading chromosome coverage matrix...")
    df = pd.read_csv(input_file, index_col=0)
    
    print(f"Original matrix shape: {df.shape}")
    print(f"Original columns: {len(df.columns)} columns")
    
    # Separate data columns from summary columns
    summary_cols = [col for col in df.columns if col in ['Mean', 'Std', 'Min', 'Max', 'Mean_Across_Samples', 'Std_Across_Samples', 'Min_Across_Samples', 'Max_Across_Samples']]
    data_cols = [col for col in df.columns if col not in summary_cols]
    
    print(f"Data columns: {len(data_cols)}")
    print(f"Summary columns: {len(summary_cols)}")
    
    # Parse and organize sample information
    samples_info = {}
    
    for col in data_cols:
        # Extract sample info: 281123_SI_14_S14_L004
        match = re.match(r'281123_([SL]I)_(\d+)_S\d+_L(\d+)', col)
        if match:
            sample_type = match.group(1)  # SI or LI
            sample_num = int(match.group(2))  # sample number
            lane_num = int(match.group(3))  # lane number
            
            sample_key = f"{sample_type}_{sample_num:02d}"  # SI_01, SI_02, etc.
            
            if sample_key not in samples_info:
                samples_info[sample_key] = {}
            
            samples_info[sample_key][lane_num] = col
    
    print(f"\nFound samples:")
    for sample_key in sorted(samples_info.keys()):
        lanes = sorted(samples_info[sample_key].keys())
        print(f"  {sample_key}: lanes {lanes}")
    
    # Create desired column order
    print("\nCreating desired column order...")
    ordered_columns = []
    
    # First SI samples (SI_01 to SI_29)
    for i in range(1, 30):
        sample_key = f"SI_{i:02d}"
        if sample_key in samples_info:
            # Add lanes in order L001, L002, L003, L004
            for lane in sorted(samples_info[sample_key].keys()):
                col_name = samples_info[sample_key][lane]
                ordered_columns.append(col_name)
                print(f"  Added: {col_name}")
    
    # Then LI samples (LI_01 to LI_04)
    for i in range(1, 5):
        sample_key = f"LI_{i:02d}"
        if sample_key in samples_info:
            # Add lanes in order L001, L002, L003, L004
            for lane in sorted(samples_info[sample_key].keys()):
                col_name = samples_info[sample_key][lane]
                ordered_columns.append(col_name)
                print(f"  Added: {col_name}")
    
    # Add any remaining columns that didn't match the pattern
    remaining_cols = [col for col in data_cols if col not in ordered_columns]
    if remaining_cols:
        print(f"\nAdding {len(remaining_cols)} remaining columns:")
        for col in remaining_cols:
            print(f"  {col}")
        ordered_columns.extend(remaining_cols)
    
    # Add summary columns at the end
    ordered_columns.extend(summary_cols)
    
    # Reorder the dataframe
    print(f"\nReordering matrix...")
    reordered_df = df[ordered_columns]
    
    print(f"New column order (first 10):")
    for i, col in enumerate(ordered_columns[:10]):
        print(f"  {i+1:2d}. {col}")
    if len(ordered_columns) > 10:
        print(f"  ... and {len(ordered_columns)-10} more")
    
    # Save reordered matrix
    reordered_df.to_csv(output_file)
    print(f"\nSaved reordered matrix to: {output_file}")
    
    # Create column order reference
    order_file = output_file.replace('.csv', '_column_order.txt')
    with open(order_file, 'w') as f:
        f.write("Column Order in Reorganized Matrix\n")
        f.write("=" * 40 + "\n\n")
        
        current_sample = None
        for i, col in enumerate(ordered_columns):
            # Extract sample info for grouping
            match = re.match(r'281123_([SL]I_\d+)', col)
            if match:
                sample = match.group(1)
                if sample != current_sample:
                    if current_sample is not None:
                        f.write("\n")
                    f.write(f"{sample}:\n")
                    current_sample = sample
                f.write(f"  {i+1:3d}. {col}\n")
            else:
                if current_sample is not None:
                    f.write("\n")
                f.write(f"Summary columns:\n")
                f.write(f"  {i+1:3d}. {col}\n")
                current_sample = None
    
    print(f"Saved column order reference to: {order_file}")
    
    # Print summary
    print(f"\n=== REORDER SUMMARY ===")
    print(f"Matrix shape: {reordered_df.shape}")
    print(f"Total columns: {len(ordered_columns)}")
    print(f"Data columns: {len(ordered_columns) - len(summary_cols)}")
    print(f"Summary columns: {len(summary_cols)}")
    
    # Show sample grouping
    print(f"\nSample grouping in new order:")
    current_sample = None
    sample_count = 0
    lane_count = 0
    
    for col in ordered_columns:
        match = re.match(r'281123_([SL]I_\d+)', col)
        if match:
            sample = match.group(1)
            if sample != current_sample:
                if current_sample is not None:
                    print(f"  {current_sample}: {lane_count} lanes")
                current_sample = sample
                sample_count += 1
                lane_count = 0
            lane_count += 1
    
    if current_sample is not None:
        print(f"  {current_sample}: {lane_count} lanes")
    
    print(f"Total samples: {sample_count}")
    
    return reordered_df

# Main execution
if __name__ == "__main__":
    # Set file paths
    input_file = "chromosome_coverage_matrix.csv"
    output_file = "chromosome_coverage_reordered.csv"
    
    print("=== CHROMOSOME MATRIX COLUMN REORDERING ===")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print()
    
    try:
        reordered_df = reorder_chromosome_matrix(input_file, output_file)
        print("\n=== SUCCESS! ===")
        print("Files created:")
        print(f"  - {output_file}")
        print(f"  - {output_file.replace('.csv', '_column_order.txt')}")
        
    except FileNotFoundError:
        print(f"ERROR: Could not find {input_file}")
        print("Make sure you're in the directory with the chromosome coverage matrix.")
    except Exception as e:
        print(f"ERROR: {e}")