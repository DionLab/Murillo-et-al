import pandas as pd
import os
import gzip
from pathlib import Path
import sys
from multiprocessing import Pool
import time

def clean_sample_name(filename):
    cleaned = filename.replace('.per-base.bed.gz', '')
    cleaned = cleaned.replace('.regions.bed.gz', '')
    cleaned = cleaned.replace('.bed.gz', '')
    cleaned = cleaned.split('.')[0]
    return cleaned

def process_single_file(file_path):
    """Process a single BED file efficiently"""
    try:
        sample_name = clean_sample_name(os.path.basename(file_path))
        print(f"  Processing {sample_name}...")
        
        chr_totals = {}
        chr_lengths = {}
        lines_processed = 0
        
        if file_path.endswith('.gz'):
            file_handle = gzip.open(file_path, 'rt')
        else:
            file_handle = open(file_path, 'r')
        
        try:
            for line in file_handle:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                try:
                    chr_name = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    coverage = float(parts[3])
                    
                    length = end - start
                    weighted_cov = coverage * length
                    
                    if chr_name not in chr_totals:
                        chr_totals[chr_name] = 0
                        chr_lengths[chr_name] = 0
                    
                    chr_totals[chr_name] += weighted_cov
                    chr_lengths[chr_name] += length
                    lines_processed += 1
                    
                except (ValueError, IndexError):
                    continue
        
        finally:
            file_handle.close()
        
        if not chr_totals:
            return None, sample_name
        
        # Calculate averages
        chr_averages = {}
        for chr_name in chr_totals:
            if chr_lengths[chr_name] > 0:
                chr_averages[chr_name] = round(chr_totals[chr_name] / chr_lengths[chr_name], 3)
        
        avg_coverage = sum(chr_averages.values()) / len(chr_averages)
        print(f"  {sample_name}: {avg_coverage:.2f}x avg, {len(chr_averages)} chrs, {lines_processed:,} lines")
        
        return pd.Series(chr_averages), sample_name
        
    except Exception as e:
        print(f"  ERROR {sample_name}: {e}")
        return None, clean_sample_name(os.path.basename(file_path))

def process_batch(file_batch, batch_num):
    """Process a batch of files using multiprocessing"""
    print(f"\n=== BATCH {batch_num}: Processing {len(file_batch)} files ===")
    start_time = time.time()
    
    # Use multiprocessing for this batch
    with Pool(processes=min(8, len(file_batch))) as pool:
        results = pool.map(process_single_file, file_batch)
    
    # Collect successful results
    batch_data = {}
    successful = 0
    failed = 0
    
    for chr_coverage, sample_name in results:
        if chr_coverage is not None:
            batch_data[sample_name] = chr_coverage
            successful += 1
        else:
            failed += 1
    
    elapsed = time.time() - start_time
    print(f"BATCH {batch_num} COMPLETE: {successful} success, {failed} failed in {elapsed:.1f}s")
    
    return batch_data

def main():
    folder_path = "/scratch/c.mpmrrp/WGS_projects/S4_WGS_dec23/output/mosdepth_windowed"
    
    # Get all BED files
    bed_files = []
    for pattern in ["*.regions.bed.gz", "*.per-base.bed.gz", "*.bed.gz"]:
        bed_files.extend(list(Path(folder_path).glob(pattern)))
    
    bed_files = list(set([str(f) for f in bed_files]))
    print(f"Found {len(bed_files)} BED files")
    
    if not bed_files:
        print("No files found!")
        return
    
    # Split files into batches of 8
    batch_size = 8
    file_batches = [bed_files[i:i+batch_size] for i in range(0, len(bed_files), batch_size)]
    
    print(f"Processing {len(file_batches)} batches of {batch_size} files each")
    
    # Process batches
    all_data = {}
    total_successful = 0
    total_failed = 0
    
    for batch_num, file_batch in enumerate(file_batches, 1):
        batch_data = process_batch(file_batch, batch_num)
        
        # Merge batch data
        all_data.update(batch_data)
        batch_successful = len(batch_data)
        batch_failed = len(file_batch) - batch_successful
        
        total_successful += batch_successful
        total_failed += batch_failed
        
        print(f"Running totals: {total_successful} successful, {total_failed} failed")
        
        # Small pause between batches to let system recover
        if batch_num < len(file_batches):
            print("Pausing 5 seconds between batches...")
            time.sleep(5)
    
    print(f"\n=== ALL BATCHES COMPLETE ===")
    print(f"Final results: {total_successful} successful, {total_failed} failed")
    
    if not all_data:
        print("No data to save!")
        return
    
    # Create final matrix
    print("Creating final matrix...")
    df = pd.DataFrame(all_data)
    
    # Sort chromosomes
    def chr_sort_key(chr_name):
        chr_name = str(chr_name)
        if chr_name.startswith('chr'):
            chr_num = chr_name[3:]
        else:
            chr_num = chr_name
        
        if chr_num.isdigit():
            return (0, int(chr_num))
        elif chr_num in ['X', 'Y']:
            return (1, ord(chr_num))
        else:
            return (2, chr_num)
    
    sorted_indices = sorted(df.index, key=chr_sort_key)
    df = df.reindex(sorted_indices)
    
    # Add summary stats
    df['Mean'] = df.mean(axis=1).round(3)
    df['Std'] = df.std(axis=1).round(3)
    df['Min'] = df.min(axis=1).round(3)
    df['Max'] = df.max(axis=1).round(3)
    
    # Save results
    print("Saving results...")
    
    # Main matrix
    df.to_csv("chromosome_coverage_matrix.csv")
    print(f"Saved: chromosome_coverage_matrix.csv ({df.shape})")
    
    # Sample summary
    sample_cols = [col for col in df.columns if col not in ['Mean', 'Std', 'Min', 'Max']]
    sample_summary = pd.DataFrame({
        'Sample': sample_cols,
        'Avg_Coverage': [df[col].mean() for col in sample_cols],
        'Std_Coverage': [df[col].std() for col in sample_cols]
    }).round(3)
    sample_summary.to_csv("sample_summary.csv", index=False)
    print(f"Saved: sample_summary.csv ({len(sample_summary)} samples)")
    
    # Chromosome summary
    chr_summary = pd.DataFrame({
        'Chromosome': df.index,
        'Average_Coverage': df['Mean'],
        'Coverage_Std': df['Std'],
        'Min_Sample_Coverage': df['Min'],
        'Max_Sample_Coverage': df['Max']
    })
    chr_summary.to_csv("chromosome_summary.csv", index=False)
    print(f"Saved: chromosome_summary.csv ({len(chr_summary)} chromosomes)")
    
    # Quick stats
    overall_avg = df[sample_cols].mean().mean()
    print(f"\nOverall average coverage: {overall_avg:.2f}x")
    print(f"Coverage range: {df[sample_cols].mean().min():.2f}x - {df[sample_cols].mean().max():.2f}x")
    print(f"Number of chromosomes: {len(df)}")
    print(f"Number of samples: {len(sample_cols)}")
    
    print("\n=== ANALYSIS COMPLETE ===")

if __name__ == "__main__":
    main()
