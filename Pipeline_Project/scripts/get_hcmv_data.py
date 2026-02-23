#this script fetches the HCMV genome from NCBI, extracts the CDS features, and saves them in a FASTA file for Snakemake to use.
from Bio import Entrez, SeqIO
import os

# this will be used by NCBI to contact us if there are any issues with our requests, so it's important to set it to a valid email address.
Entrez.email = "pfay@luc.edu" 
#here i set the API key to increase the rate limit for NCBI requests. You should replace "your_api_key_here" with your actual NCBI API key if you have one. 
# If you don't have an API key, you can leave this line out
def fetch_data():
    # NC_006273.2 is the exact RefSeq ID for the HCMV genome
    accession = "NC_006273.2" 
    #this is the main function that fetches the data, extracts the CDS features, and saves them in a FASTA file. 
    # It uses Biopython's Entrez module to fetch the GenBank record and SeqIO to parse it.
    print(f"Downloading {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    #i did this to extract the CDS features from the GenBank record. For each CDS feature, we get the protein_id qualifier to use as the FASTA header and extract the corresponding sequence. 
    # We then create a new SeqRecord for each CDS and store them in a list.
    cds_records = []
    for feature in record.features:
        if feature.type == "CDS":
            # We need the protein_id for the fasta header
            prot_id = feature.qualifiers.get("protein_id", ["no_id"])[0]
            sequence = feature.extract(record.seq)
            
            from Bio.SeqRecord import SeqRecord
            new_record = SeqRecord(sequence, id=prot_id, description="")
            cds_records.append(new_record)

    # Save the file Snakemake is looking for
    #this will write the extracted CDS sequences to a FASTA file named "hcmv_cds.fasta". Each sequence will have its protein_id as the header.
    SeqIO.write(cds_records, "hcmv_cds.fasta", "fasta")
    print(f"Success: Extracted {len(cds_records)} CDS features.")

# lastly, this block checks if the script is being run directly (as the main program) and calls the fetch_data function to execute the data fetching and processing steps.
if __name__ == "__main__":
    fetch_data()