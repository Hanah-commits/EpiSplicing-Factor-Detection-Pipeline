import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse
from collections import Counter
import json


def is_valid_url(url):
    # Check if the URL is valid using urlparse
    parsed = urlparse(url)
    return all([parsed.scheme, parsed.netloc])


def convert_to_iupac(motifs):

    iupac_codes = {
        'R':['A', 'G'],
        'Y':['C','U'],
        'S':['G','C'],
        'W':['A','U'],
        'K':['G','U'],
        'M':['A','C']
    }

    output_seq = []
    
    for base1, base2 in zip(motifs[0], motifs[1]):
        if base1 == base2:
            output_seq.append(base1)
        else:
            # get IUPAC code
            for code, bases in iupac_codes.items():
                if sorted([base1, base2]) == sorted(bases):
                    output_seq.append(code)
                    break

    return ''.join(output_seq)


def consensus_sequence(sequences):

    transposed = zip(*sequences)
    
    consensus = []
    for position in transposed: # iterate over positions
        # count the frequency of each nucleotide at the current position
        frequency = Counter(position)
        
        # get the nucleotide with the highest frequency
        consensus_nucleotide = frequency.most_common(1)[0][0]
        consensus.append(consensus_nucleotide)
    
    return ''.join(consensus)


def crawl_db():
    rbps = ["AGGF1", "AKAP1", "AQR", "BUD13", "CSTF2T", "DDX3X", "DDX52", "DDX55", "DDX6", "DGCR8", "DHX30", "DROSHA", "EFTUD2", "EXOSC5", "FAM120A", "FASTKD2", "FTO", "FUS", "FXR2", "GRWD1", "GTF2F1", "HLTF", "HNRNPA1", "HNRNPC", "HNRNPK", "HNRNPL", "HNRNPM", "HNRNPU", "HNRNPUL1", "IGF2BP1", "ILF3", "KHSRP", "LARP4", "LARP7", "LIN28B", "LSM11", "MATR3", "NCBP2", "NOLC1", "NSUN2", "PCBP1", "PRPF8", "PTBP1", "QKI", "RBFOX2", "RBM15", "RBM22", "RPS3", "SAFB", "SDAD1", "SF3B4", "SLTM", "SMNDC1", "SND1", "SRSF1", "SRSF7", "SRSF9", "SSB", "SUPV3L1", "TAF15", "TARDBP", "TBRG4", "TIA1", "TRA2A", "TROVE2", "U2AF1", "U2AF2", "UCHL5", "UPF1", "UTP18", "WDR43", "XRCC6", "XRN2", "YBX3", "ZC3H11A", "ZNF800"]
    all_rbp_motifs = {}

    for rbp in rbps:
        rbp_motifs = []
        for cell in ['HepG2','K562']:
            url= f'https://zhanglab.c2b2.columbia.edu/mCrossBase/rbp.php?id={cell}.{rbp}'
            try:
                response = requests.get(url)
                response.raise_for_status()  # exception for HTTP errors
            except requests.exceptions.RequestException as e:
                print(f"Failed to retrieve the URL: {e} for {rbp}")
                continue

           # parse content
            soup = BeautifulSoup(response.text, 'html.parser')
            
            print(f'Processing: {rbp} - {cell}')
            
            try:
                rows = soup.find_all('tr', {'class': 'align-text-bottom'})
            except:
                print(f"No motifs available: {e} for {rbp}")
                continue

            #get consensus motifs
            for row in rows:
                cell = row.find_all('td')[2]
                num_sites = int(row.find_all('td')[5].text)
                if num_sites > 100: # min 50 binding sites supporting this motif
                    motifs = []
                    for li in cell.find_all("li"): 
                        motifs.append(li.text.replace("T", "U"))
                    
                    if len(motifs) ==2:
                        motifs = convert_to_iupac(motifs)
                        rbp_motifs.append(motifs)
                    elif len(motifs) > 2:
                        motifs = consensus_sequence(motifs)
                        rbp_motifs.append(motifs)
                    else:
                        rbp_motifs.append(motifs[0])
                else:
                    'Insufficient binding sites to support motif'

        #drop duplicates
        all_rbp_motifs[rbp] = list(set(rbp_motifs))
            
    with open('rbp_motifs.json', 'w') as fp:
        json.dump(all_rbp_motifs, fp, indent=4)


if __name__ == '__main__':
    crawl_db()




