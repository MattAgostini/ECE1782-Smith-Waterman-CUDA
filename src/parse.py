import os
import math

datafile = open("uniprot_subset.dat", 'r')

'''
SQ   SEQUENCE   256 AA;  29735 MW;  B4840739BF7D4121 CRC64;
     MAFSAEDVLK EYDRRRRMEA LLLSLYYPND RKLLDYKEWS PPRVQVECPK APVEWNNPPS
     EKGLIVGHFS GIKYKGEKAQ ASEVDVNKMC CWVSKFKDAM RRYQGIQTCK IPGKVLSDLD
     AKIKAYNLTV EGVEGFVRYS RVTKQHVAAF LKELRHSKQY ENVNLIHYIL TDKRVDIQHL
     EKDLVKDFKA LVESAHRMRQ GHMINVKYIL YQLLKKHGHG PDGPDILTVK TGSKGVLYDD
     SFRKIYTDLG WKFTPL
'''

# Parameter obtained from the format
sequenceWidth = 60.0
             
sequence = []

line = datafile.readline()

# Find sequences and associated sequence lengths
while line:
    if (line.startswith("SQ")):
        seqLen = int(line.split()[2])

        parsedSequence = ""
        lines = math.ceil(seqLen / sequenceWidth)
        for i in range(int(lines)):
            unparsedSequence = (datafile.readline()).split()
            for j in range(len(unparsedSequence)):
                parsedSequence += unparsedSequence[j]
                
        #print parsedSequence + "\n"
        sequence.append(parsedSequence)

    line = datafile.readline()

# Sort sequence by length
sequence.sort(key = lambda s: len(s))

parsedfile = open("uniprot_subset_p.dat", 'w+')

for i in range(len(sequence)):
    #print sequence[i]
    parsedfile.write(sequence[i] + "\n")








