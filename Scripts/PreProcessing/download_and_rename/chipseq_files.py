import pandas as pd
import os
import sys
import requests
import subprocess


def create_url(file_accession, file_format):
    base_url = 'https://www.encodeproject.org/files/'
    if file_format == 'bam':
        return f'{base_url}{file_accession}/@@download/{file_accession}.bam'
    elif file_format == 'bed narrowPeak':
        return f'{base_url}{file_accession}/@@download/{file_accession}.bed.gz'


def download_file(url, local_filename):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
                print(f"Downloading {local_filename}: {f.tell() / total_size:.2%}", end='\r')
    print(f"Downloaded {local_filename}")
    return local_filename


def merge_bam_files(file_accessions, output_file):
    input_files = " ".join([f"{fa}.bam" for fa in file_accessions])
    cmd = f'samtools merge -f {output_file} {input_files}'
    print(f"Executing: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"Merged BAM files into {output_file}")

def bam_to_bed(input_file, output_file):
    cmd = f'bedtools bamtobed -i {input_file} > {output_file}'
    print(f"Executing: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"Converted {input_file} to {output_file}")


def sort_and_merge_bed(input_file):
    cmd = f'sort -k1,1 -k2,2n {input_file} | bedtools merge -i - > tmp.bed && mv tmp.bed {input_file}'
    print(f"Executing: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    print(f"Sorted and merged peaks in {input_file}")


if __name__ == "__main__":
    # Load the data
    df = pd.read_csv(sys.argv[1], sep='\t')

    # Group by 'Biosample', 'File format', and 'Experiment target'
    grouped = df.groupby(['Biosample', df['File format'] == 'bam', 'Experiment target'])

    for (biosample, is_bam, target), group in grouped:
        if is_bam:
            # Process BAM files

            # Step 1: Prepare link
            file_accessions = group['File accession'].tolist()
            urls = [create_url(fa, 'bam') for fa in file_accessions]
            
            # Step 2: Download BAM files
            for url in urls:
                download_file(url, url.split('/')[-1])
            
            # Step 3: Merge BAM files
            output_bam = f'{target}_{biosample}_alignment.bam'
            merge_bam_files(file_accessions, output_bam)
            
            # Step 4: Convert to BED format
            output_bed = f'{output_bam.replace(".bam", ".bed")}'
            bam_to_bed(output_bam, output_bed)
            
            # # Delete BAM files
            # for fa in file_accessions:
            #     os.remove(f'{fa}.bam')
            #     print(f"Deleted {fa}.bam")

        else:
            # Process narrowPeak files

            # Step 1: Prepare link
            file_accessions = group['File accession'].tolist()
            urls = [create_url(fa, 'bed narrowPeak') for fa in file_accessions]
            
            # Step 2: Download and unzip BED files
            for url in urls:
                local_file = download_file(url, url.split('/')[-1])
                cmd = f'gzip -d {local_file}'
                print(f"Executing: {cmd}")
                subprocess.run(cmd, shell=True, check=True)
                unzipped_file = local_file.replace('.gz', '')
                print(f"Unzipped {local_file} to {unzipped_file}")
                
                # Step 3: Sort and merge peaks
                sort_and_merge_bed(unzipped_file)
                
                # Step 4: Rename file
                final_bed = f'{target}_{biosample}_peak.bed'
                os.rename(unzipped_file, final_bed)
                print(f"Renamed {unzipped_file} to {final_bed}")
