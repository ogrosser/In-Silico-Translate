from random import choice
import re #Allows for use of regular expressions

def dna(length): #Function that generates a random nucleotide sequence of specified length
    return ''.join(choice("ATGC") for i in range(length))

ten_seq = []
for i in range(10): #Create list of 10 nucleotide sequences
    ten_seq.append(dna(200))

def r_comp(seq): #Returns reverse complement of nucleotide sequence
    nuc = {'A':'T', 'T':'A', 'C':'G', 'G':'C'} #nucleotide pairs
    return ''.join(nuc[n] for n in reversed(seq))

o = open('protein.txt', 'w') #Create output text file. Overwrites file if it already exists
o.close()

#Dictionary of all codons and their corresponding protein codes
proteins = {'ATG':'M', 'ATA': 'I', 'ATC':'I', 'ATT':'I',
            'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
            'AAT': 'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
            'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
            'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
            'TAT':'Y', 'TAC':'Y', 'TAA':'-', 'TAG':'-',
            'TGT':'C', 'TGC':'C', 'TGA':'-', 'TGG':'W',
            'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
            'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R'}

def find_codons(seq): #Find indexes of all instances of start and stop codons
    global start_index
    start_index = []
    global stop_index
    stop_index = {}
    start = re.finditer(r"ATG", seq) #All indexes of start codon
    stop = re.finditer(r"TAG|TGA|TAA", seq) #All indexes of stop codons
    if start is None or stop is None: #If there are either no start or stop codons
        f.write("ORF not found due to missing start or stop codon.\n")
    else:
        for m in start: #start_index contains all indexes of start codon ATG
            start_index.append(m.start())
        for n in stop: #stop_index contains all indexes of stop codons TAG, TGA, and TAA
        ##If codon key is already in dictionary, add index to value list. If not, then create value list
        ##with current index
            if n.group() in stop_index:
                stop_index[n.group()].append(n.start())
            else:
                stop_index[n.group()] = [n.start()]

def indiv_sil(seq): #Returns corresponding protein sequence of a single nucleotide sequence
    pro_codons = []
    pro_seq = []
    y = []
    #For each instance of ATG, list "y" contains all indexes of stop codons that fit alongside the
    #start codon when following with no overlaps against previous codons
    for i in start_index:
        y = [] #List of stop codons that fit in reading frame with current ATG
        for key in stop_index:
            for val in stop_index[key]:
                #If the distance between start and stop codon can be divided by 3 with no remainder
                if (val - i) % 3 == 0 and val > i:
                    y.append(val)
        if not y: #If current ATG has no stop codons downstream that fit, move to next ATG
            continue
        else:
            y.sort()
            #Nucleotide sequence, starting from start codon and ending before closest stop codon, split into codons 
            pro_codons = re.finditer(r'...', seq[i:y[0]])
            #Proteins corresponding to above codons. This is a nested list with each protein sequence
            #in separate element 
            pro_seq.append([proteins[p.group()] for p in pro_codons])
    if not pro_seq: #If codons do not work with current start codon resulting in no translated proteins
        f.write("No ORF found. No stop codons fit in with any start codons.\n")
    else:
        f.write("Protein sequences:\n")
        for i in range(len(pro_seq)): #Join protein codes in each element together
            f.write(''.join(str(pro_seq[i][j]) for j in range(len(pro_seq[i]))) + "\n")

def in_silico(seq_list=[]): #in-silico translation of nucleotide sequences and their reverse complement
    global f #F variable is global to allow other functions to write to output file
    with open('protein.txt', 'a') as f:
        for s in seq_list:
            f.write("\nNucleotide sequence:\n" + s + "\n")
            find_codons(s)
            indiv_sil(s)
            f.write("Reverse Complement Translation\n")
            rc = r_comp(s)
            find_codons(rc)
            indiv_sil(rc)
    f.close()

in_silico(ten_seq)
