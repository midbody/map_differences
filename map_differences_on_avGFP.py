#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 22:53:54 2024

@author: nori
"""

import os
from Bio import SeqIO, AlignIO
import shutil
from pathlib import Path
import subprocess

#### 1. Set the sequences to compare. Use the names on FPbase ####
"""
name1 = "meffGFP"
name2 = "meffRFP"
diff_color = "orange"

name1 = "rfloGFP"
name2 = "rfloRFP"
diff_color = "darkturquoise"

name1 = "ZsGreen"
name2 = "zoan2RFP"
diff_color = "salmon"

name1 = "ccalGFP1"
name2 = "ccalRFP1"
diff_color = "darkseagreen"

name1 = "mkG"
name2 = "KO"
diff_color = "palevioletred"

name1 = "cmFP512"
name2 = "OFP"
diff_color = "tan"
"""

name1 = "rfloGFP"
name2 = "rfloRFP"
diff_color = "darkturquoise"




#### 2. Set the fixed parameters
root_dir_path = "/Users/nori/Dropbox/BBSRC_EvoFPs/FP_phylogenetics/map_differences"
all_natural_proteins_fasta_path = "/Users/nori/Dropbox/BBSRC_EvoFPs/FP_phylogenetics/FPs_all_natural/natural_proteins_sequences.fasta"
avGFP_pdb_path = "/Users/nori/Dropbox/StayRed/PDB_files/1ema_posed.pdb"
name_ref = "avGFP"


#### 3. Prepare the directory to save the results ####
target = f"{name1}_vs_{name2}"
root_dir = Path(root_dir_path)
target_dir = root_dir / target
png_path = target_dir / f"{target}.png"
out_script_path = target_dir / f"{name1}_vs_{name2}.cxc"


# Check if the folder exists
if os.path.exists(target_dir):
    # If it exists, change to that directory
    os.chdir(target_dir)
    print(f"Changed to directory: {target_dir}")
else:
    # If it doesn't exist, create the folder and then change to it
    os.makedirs(target_dir)
    os.chdir(target_dir)
    print(f"Created and changed to directory: {target_dir}")


#### 4. Retrieve the sequences  ####

def get_sequence_by_protein_name(fasta_file_path, target_name):
    # Loop through each sequence record in the fasta file
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        # Split the header into parts using "|" and "/"
        header_parts = record.description.split("|")
        
        # Extract protein name (assumed to be the second part of the header)
        protein_name = header_parts[1]
        
        # Check if the protein name matches the target name
        if protein_name == target_name:
            # Return the sequence as a string if there is a match
            return str(record.seq)
    
    # Return None if the protein name was not found
    return None

# read the sequences by name
seq1 = get_sequence_by_protein_name(all_natural_proteins_fasta_path, name1)
seq2 = get_sequence_by_protein_name(all_natural_proteins_fasta_path, name2)
seq_ref = get_sequence_by_protein_name(all_natural_proteins_fasta_path, name_ref)


# Path to input and output files
input_fasta_file = f"{target}.fasta"
aligned_fasta_file = f"{target}_aligned.fasta"

#### 5. Save the sequences in a FAST file, let  MAFFT to make the aligned sequences ####

# Write sequences to a FASTA file for alignment
with open(input_fasta_file, "w") as f:
    f.write(f">{name1}\n" + seq1 + "\n")
    f.write(f">{name2}\n" + seq2 + "\n")
    f.write(">{name_ref}\n" + seq_ref + "\n")

mafft_path = shutil.which("mafft")
bash_command = f"{mafft_path} --auto {input_fasta_file} > {aligned_fasta_file}"
subprocess.run(bash_command, shell=True)


#### 6. Read the algined sequences and find the residues different between sequence 1 and sequence 2  ####

# Read the aligned sequences
alignment = AlignIO.read(aligned_fasta_file, "fasta")

# Ensure sequences A, B, and C are in the correct order in the alignment
seq1_aligned = alignment[0].seq
seq2_aligned = alignment[1].seq
ref_aligned = alignment[2].seq

# Find positions where A and B differ and use residue positions from C
differences = []
for i, (res1, res2, res_ref) in enumerate(zip(seq1_aligned, seq2_aligned, ref_aligned)):
    if res1 != res2 and res_ref != "-":  # Ignore gaps
        # Convert alignment index to residue number in C
        residue_number_in_ref = sum(1 for j in range(i + 1) if ref_aligned[j] != "-")
        differences.append((residue_number_in_ref, res1, res2))

different_residues = [item[0] for item in differences]

residues_string = ",".join(map(str, different_residues))



with open(out_script_path, "w") as script_file:
    ## set the window setting
    script_file.write("close session\n")
    script_file.write("set bgColor white\n")
    script_file.write("windowsize 1000 1000\n")
    
    ## open 1EMA (avGFP) and adjust the representation
    script_file.write(f"open {avGFP_pdb_path}\n")
    script_file.write("zoom 1.2\n")
    script_file.write("turn y -30\n")
    script_file.write("color #1 gray\n")
    script_file.write("display #1:65-67\n")
    script_file.write("~ribbon #1:64-68\n")
    script_file.write("display #1:64@N,CA,C,O\n")
    script_file.write("display #1:68@N,CA,C,O\n")  
    script_file.write("transparency #1 80 ribbons\n")
    

    ## color residues different between seq1 and seq2  
    script_file.write(f"color #1:{residues_string} {diff_color}\n")
    script_file.write(f"display #1:{residues_string}\n")
    script_file.write(f"transparency #1:{residues_string} 20 ribbons\n")
    
    ## check the mutation in the chromophore triplet 
    script_file.write("color #1:64-68 lightgray\n")
    if 65 in different_residues:
        script_file.write(f"color #1:66@N1,CA1,CB1,CG1,OG1,C1 {diff_color}\n")
        
    script_file.write("color #1 byhet\n")

    ## save a png image
    script_file.write(f"save {png_path} supersample 3")
    
   
print(f"{name1} vs {name2}")