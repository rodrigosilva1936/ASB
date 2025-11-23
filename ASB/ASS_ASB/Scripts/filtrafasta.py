from Bio import SeqIO
import re


def clean_headers(input_file, output_file):
    """
    Remove accession numbers from FASTA headers, keeping only 'Anopheles...' part
    """
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            # Remove accession number (everything before 'Anopheles')
            new_id = re.sub(r'^.*?(Anopheles.*)', r'\1', record.description)
            
            # Update record ID and description
            record.id = new_id.split()[0]  # Use first word as ID
            record.description = new_id
            
            # Write modified record
            SeqIO.write(record, out_handle, "fasta")

if __name__ == "__main__":
    
    print("Usage: python clean_headers.py input.fasta output.fasta")
       
    
    input_fasta = "C:/Users/rgs12/Downloads/align_white2.fasta"
    output_fasta = "C:/Users/rgs12/Downloads/filtered_white.fasta"
    
    clean_headers(input_fasta, output_fasta)
    print(f"Cleaned headers saved to {output_fasta}")
