import pandas as pd
import requests, sys
import subprocess
from pathlib import Path



def create_url(file_accession):
    return f'https://www.encodeproject.org/files/{file_accession}/@@download/{file_accession}.bed.gz'


def download_file(url, local_filename):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
                print(f"Downloading {local_filename}: {f.tell() / int(r.headers.get('content-length', 1)):.2%}", end='\r')
    print(f"\nDownloaded {local_filename}")
    return local_filename


if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], sep='\t')
    df = df[df.BioSample.notna()]

    # Iterate over the DataFrame and process each file
    for index, row in df.iterrows():
        bio_sample = row['BioSample']
        experiment_target = row['Experiment target']
        file_accession = row['File accession']
        
        # Step 1: Prepare link
        url = create_url(file_accession)
        
        # Step 2: Download file
        local_filename = download_file(url, f"{file_accession}.bed.gz")
        
        # Step 3: Unzip the file
        subprocess.run(f'gzip -d {local_filename}', shell=True, check=True)
        unzipped_file = local_filename.replace('.gz', '')
        
        # Step 4: Rename the file
        new_filename = f"{experiment_target}_{bio_sample}.bed"
        Path(unzipped_file).rename(new_filename)
        print(f"Renamed {unzipped_file} to {new_filename}")
