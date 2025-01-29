# In-Silico-Translate
This project is a from-scratch DNA to protein translator made without using Biopython. The translation will skip transcribing the DNA sequnces to mRNA, and instead translate DNA directly to protein.

To execute script, type the following in your terminal:
    python3 in_silico.py

The script will write each random DNA sequence and its protein sequence to this text file:
    protein.txt
    
The script will create 10 random DNA sequences 200bp long, and translate any ORFs (including nested ORFs) to their corresponding protein sequences. The translation will occur at both 5'-3' and 3'-5'. If no ORF is found, in_silico.py will write "ORF not found" and the specific reason why to protein.txt.
