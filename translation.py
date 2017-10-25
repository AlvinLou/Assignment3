import re
DNA1 = []
ORF = []
find = []
i = 0
dna2aa = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "", "TAG": "", "TGT": "C", "TGC": "C", "TGA": "", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
def translation(seq1):
    l = seq1
    proteinsequence = ""
    for i in range(0, l - 1, 3):
        coden = ORF[i:i + 3]
        if dna2aa[coden] == '':
            return proteinsequence
        proteinsequence += dna2aa[coden]
    return proteinsequence
with open("./DNA1.txt", 'r') as f:
    for line in f:
        DNA1.append(line)
        for i in range(len(DNA1)):
            a = re.search('(ATG[ATCG]*(TAG|TAA|TGA))',DNA1[i])
            if a:
                ORF = a.group(1)
                print("The DNA sequence is: "+a.group(i))
                b = len(ORF)
                c = b//3
                d = b - 3 * c
                if d == 0 :
                    protein1 = translation(len(a.group(i)))
                    print("The protein sequence is: "+protein1+"\n")
                else:
                    break
            g = str(a.group(i))
        for line in g:
            for i in range(len(DNA1)):
                i = 1
                a = re.search('(ATG[ATCG]*(TAG|TAA|TGA))', g)
                if a:
                    ORF = a.group(i)
                    b = len(ORF)
                    c = b // 3
                    d = b - 3 * c
                    if d == 0:
                        protein1 = translation(len(a.group(i)))
                        print("The DNA sequence is: " + a.group(i))
                        print("The protein sequence is: " + protein1+"\n")
                        i = i + 1
                    else:
                        break
            g = str(g)[i:]