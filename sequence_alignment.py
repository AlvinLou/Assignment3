import random

dna2aa = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"stop","TAG":"stop","TGT":"C","TGC":"C","TGA":"stop","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}
aai = {
    "A":0,"R":1,"N":2,"D":3,"C":4,"Q":5,"E":6,"G":7,"H":8,"I":9,
    "L":10,"K":11,"M":12,"F":13,"P":14,"S":15,"T":16,"W":17,"Y":18,"V":19
}
blosum = [[0]*20 for i in range(20)]
blosum[0] = [4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[1] = [-1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[2] = [-2, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[3] = [-2, -2, 1, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[4] = [0, -3, -3, -3, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[5] = [-1, 1, 0, 0, -3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[6] = [-1, 0, 0, 2, -4, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[7] = [0, -2, 0, -1, -3, -2, -2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[8] = [-2, 0, 1, -1, -3, 0, 0, -2, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[9] = [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[10] = [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[11] = [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, 0, 0, 0, 0, 0, 0, 0, 0]
blosum[12] = [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, 0, 0, 0, 0, 0, 0]
blosum[13] = [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, 0, 0, 0, 0, 0, 0]
blosum[14] = [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, 0, 0, 0, 0, 0]
blosum[15] = [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 0, 0, 0, 0]
blosum[16] = [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, 0, 0, 0]
blosum[17] = [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 0, 0]
blosum[18] = [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, 0]
blosum[19] = [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4]
for i in range(20):
    for j in range(i+1,20):
        blosum[i][j] = blosum[j][i]
#for i in range(20):
    #print(blosum[i])

def translation(seq1):
    l = len(seq1)
    pro = ""
    for i in range(0,l-1,3):
        coden = seq1[i:i+3]
        if dna2aa[coden] == '':
            return pro
        if dna2aa[coden] == 'stop':
            break
        pro += dna2aa[coden]
        #print("check point1:"+coden)
    return pro

def Protein_local_align(seq1, seq2):
    print("The length of sequence1:%d" % len(seq1))
    print("The length of sequence2:%d" % len(seq2))
    l1 = len(seq1) + 1
    l2 = len(seq2) + 1
    gap = -2
    #mismatch = -2
    #match = 3
    m = [[0] * l1 for i in range(l2)]
    d = [['X'] * l1 for i in range(l2)]
    for i in range(1, l1):
        m[0][i] = 0  # m[0][i-1] + gap
        # d[0][i] = "L"

    for j in range(1, l2):
        m[j][0] = 0  # m[j-1][0] + gap
        # d[j][0] = "U"

    max_number = -1000
    x = 0
    y = 0
    for i in range(1, l2):
        for j in range(1, l1):
            pt_u = m[i - 1][j] + gap
            pt_l = m[i][j - 1] + gap
            pt_d = m[i - 1][j - 1] + blosum[aai[seq2[i-1]]][aai[seq1[j-1]]]
            #if seq1[j - 1] == seq2[i - 1]:
            #    pt_d = m[i - 1][j - 1] + match
            #else:
            #    pt_d = m[i - 1][j - 1] + mismatch
            m[i][j] = max(pt_u, pt_l, pt_d, 0)
            if m[i][j] == pt_d:
                d[i][j] = 'D'
            elif m[i][j] == pt_u:
                d[i][j] = 'U'
            elif m[i][j] == pt_l:  # else:
                d[i][j] = 'L'

            if max_number < m[i][j]:
                max_number = m[i][j]
                x = i
                y = j
    print("The biggest number: %d " % max_number + "@ %d line and " % x + " %d column." % y)
    # for j in range(l2):
    #    print(m[j])
    alg1 = ''
    alg2 = ''
    alg0 = ''

    while (m[x][y] > 0):
        if d[x][y] == 'D':
            alg1 = seq1[y - 1] + alg1
            alg2 = seq2[x - 1] + alg2
            if seq1[y - 1] == seq2[x - 1]:
                alg0 = '|' + alg0
            elif blosum[aai[seq2[x-1]]][aai[seq1[y-1]]] > 1: #seq1[y - 1] == seq2[x - 1]:
                alg0 = ':' + alg0
            elif blosum[aai[seq2[x-1]]][aai[seq1[y-1]]] > 0: #seq1[y - 1] == seq2[x - 1]:
                alg0 = '.' + alg0
            else:
                alg0 = ' ' + alg0
            x -= 1
            y -= 1
        elif d[x][y] == 'U':
            alg1 = '-' + alg1
            alg2 = seq2[x - 1] + alg2
            alg0 = ' ' + alg0
            x -= 1
        elif d[x][y] == 'L':
            alg1 = seq1[y - 1] + alg1
            alg2 = '-' + alg2
            alg0 = ' ' + alg0
            y -= 1
        else:
            x -= 1
            y -= 1
            # print ("check point 1: %d\n"%m[x][y]+ alg1 + "\n" + alg0 + "\n" + alg2 + "x=%d "%x+"y=%d "%y)
    print("local alignment:\n" + alg1 + "\n" + alg0 + "\n" + alg2 + "\nAnd the score = %d " % max_number)


def DNA_local_align(seq1, seq2):
    print("The length of sequence1:%d" % len(seq1))
    print("The length of sequence2:%d" % len(seq2))
    l1 = len(seq1) + 1
    l2 = len(seq2) + 1
    gap = -1
    mismatch = -2
    match = 3
    m = [[0] * l1 for i in range(l2)]
    d = [['X'] * l1 for i in range(l2)]
    for i in range(1, l1):
        m[0][i] = 0  # m[0][i-1] + gap
        # d[0][i] = "L"

    for j in range(1, l2):
        m[j][0] = 0  # m[j-1][0] + gap
        # d[j][0] = "U"

    max_number = -1000
    x = 0
    y = 0
    for i in range(1, l2):
        for j in range(1, l1):
            pt_u = m[i - 1][j] + gap
            pt_l = m[i][j - 1] + gap
            if seq1[j - 1] == seq2[i - 1]:
                pt_d = m[i - 1][j - 1] + match
            else:
                pt_d = m[i - 1][j - 1] + mismatch
            m[i][j] = max(pt_u, pt_l, pt_d, 0)
            if m[i][j] == pt_d:
                d[i][j] = 'D'
            elif m[i][j] == pt_u:
                d[i][j] = 'U'
            elif m[i][j] == pt_l:  # else:
                d[i][j] = 'L'

            if max_number < m[i][j]:
                max_number = m[i][j]
                x = i
                y = j
    print("The biggest number: %d " % max_number + "@ %d line and " % x + " %d column." % y)
    # for j in range(l2):
    #    print(m[j])
    alg1 = ''
    alg2 = ''
    alg0 = ''

    while (m[x][y] > 0):
        if d[x][y] == 'D':
            alg1 = seq1[y - 1] + alg1
            alg2 = seq2[x - 1] + alg2
            if seq1[y - 1] == seq2[x - 1]:
                alg0 = '|' + alg0
            else:
                alg0 = ' ' + alg0
            x -= 1
            y -= 1
        elif d[x][y] == 'U':
            alg1 = '-' + alg1
            alg2 = seq2[x - 1] + alg2
            alg0 = ' ' + alg0
            x -= 1
        elif d[x][y] == 'L':
            alg1 = seq1[y - 1] + alg1
            alg2 = '-' + alg2
            alg0 = ' ' + alg0
            y -= 1
        else:
            x -= 1
            y -= 1
            # print ("check point 1: %d\n"%m[x][y]+ alg1 + "\n" + alg0 + "\n" + alg2 + "x=%d "%x+"y=%d "%y)
    print("local alignment:\n" + alg1 + "\n" + alg0 + "\n" + alg2 + "\nAnd the score = %d " % max_number)

def main():
    dna1 = "GCTAACTGAGCTAACTGAGCGAGCACTGGCCTGGCCTGCTGCTGCCTGTGCCATGGCTCCTGGGAGTGTCTCCAGTGTTTCTTCCTCCTCTTTTCCCTCCAGGGACACATCCCCTTCTGGATCATGTGGGCTCCCTGGAGCTGACAAGCCAGGTCCAAGTTGCCGCAGAATCCAAGCAGGCCAAAGGAACCCAACAATGCTGCACATGGTGCTAGAGGCTTTGAAGGCCCGGGAGGCACGCCAGGGCACATCAGTTGTAGCCATCAAGGTCTACATCCAACACAAGTACCCGACAGTGGACACCACCCGTTTCAAGTACCTGTTGAAGCAAGCTCTGGAAACTGGCGTTCGTCGAGGCCTCCTCACCAGGCCTGCTCACTCCAAGGCCAAGGGTGCCACTGGCAGCTTCAAACTAGTTCCAAAGCCCAAGACAAAGAAAGCCTGTGCCCCCAAAGCCGGCAGGGGAGCTGCAGGTGCCAAGGAGACAGGCTCCAAGAAATCTGGATTGCTGAAGAAAGACCAAGTTGGCAAGGCCACGATGGAGAAAGGGCAGAAGAGGAGGGCTTACCCTTGCAAGGCAGCCACACTGGAGATGGCACCTAAGAAAGCCAAGGCGAAACCGAAAGAGGTCAGAAAGGCTCCCCTAAAACAAGACAAAGCAGCAGGGGCCCCTCTGACTGCCAATGGAGGCCAGAAGGTCAAACGCAGTGGGAGCAGGCAAGAAGCAAATGCCCATGGGAAAACCAAAGGTGAGAAATCGAAGCCCTTGGCCAGCAAGGTCCAGAATAGCGTTGCTTCCCTCGCCAAAAGGAAGATGGCAGACATGGCCCACACTGTGACAGTTGTTCAGGGGGCTGAGACAGTACAGGAGACCAAAGTGCCCACTCCTTCCCAGGACATAGGACACAAAGTACAACCCATACCTAGGGTCAGGAAGGCAAAGACCCCTGAGAACACTCAGGCCTGAGTTACTTCCCAAGACCTCCTCCAAGGCTCCCAGCAAGAAGGCTGAGGCTAGTAGCTAGGGCCAGGGCTGGGGAGATGGCGATTCTGAAGCTTTTATTGTCTAATAAGCTGTCACAATGTTTAATCATAATTTATCAATAAAGACTTTGTATTCACACACTGGTA"
    dna2 = "ACTTATTGTCTTTTCTGGGAAGACAAAAACATGTCGGAGACTGCTCCACTTGCTCCTACCATTCCTGCACCCGCAGAAAAAACACCTGTGAAGAAAAAGGCGAAGAAGGCAGGCGCAACTGCTGGGAAACGCAAAGCATCCGGACCCCCAGTATCTGAGCTTATCACCAAGGCAGTGGCAGCTTCTAAGGAGCGCAGCGGCGTTTCTCTGGCCGCGCTTAAGAAAGCGCTTGCGGCTGCTGGCTACGATGTAGAAAAAAACAACAGCCGTATCAAGCTTGGCCTCAAGAGCTTGGTGAGCAAAGGTACTCTGGTGCAGACCAAAGGTACCGGTGCTTCTGGCTCCTTCAAACTCAACAAGAAAGCGGCTTCCGGGGAAGGCAAACCCAAGGCCAAAAAGGCTGGCGCAGCCAAGCCTAGGAAGCCTGCTGGGGCAGCCAAGAAGCCCAAGAAGGTGGCTGGCGCCGCTACCCCGAAGAAAAGCATCAAAAAGACTCCTAAGAAGGTAAAGAAGCCAGCAACCGCTGCTGGGACCAAGAAAGTGGCCAAGAGTGCGAAAAAGGTGAAAACACCTCAGCCAAAAAAAGCTGCCAAGAGTCCAGCTAAGGCCAAAGCCCCTAAGCCCAAGGCGGCCAAGCCTAAGTCGGGGAAGCCGAAGGTTACAAAGGCAAAGAAGGCAGCTCCGAAGAAAAAGTGAAACTGGCGGGACGTTCCCCTTTGAAAATTTTAAACGGCTCTTTTCAGAGCCACCCA®"

    DNA_local_align(dna1, dna2)
    for i in range(3):
        input_seq = dna1[i:]
        protein1 = translation(input_seq)
        for j in range(3):
            input_seq = dna2[j:]
            protein2 = translation(input_seq)
            #print("check point2:"+input_seq)
            print("protein 1 frame %d: " %i+protein1+"\nprotein 2 frame %d: " %j+protein2)

            Protein_local_align(protein1, protein2)
main()


def random_seq(seq1):
    l = len(seq1)
    n = 0
    m = 0
    for i in range(l):
        n = random.randint(0, l)
        m = n
        while (n == m):
            m = random.randint(0, l)
        temp = seq1[n]
        seq1[n] = seq1[m]
        seq1[m] = temp
    return seq1


def calculate_average(score):
    l = len(score)
    sum = 0
    for i in range(l):
        sum += score[i]


def calculate_std_dev(score, avg):
    l = len(score)
    variance = 0.0
    for i in range(1):
        variance += (float(score[i]) - avg) ** 2
        std_dev = (variance/1) ** .5

#def stat(seq1, seq2):
 #    no_interaction = 1000
  #   pro_align_score = protein_loacal_align(seq1, seq2)
  #   Score = []
  #   for I in range(no_interaction):
  #      Seq3 = random_seq(seq1)
  #      seq4 = random_seq(seq2)
  #      score[i] = protein_loacal_align(seq3, seq4)
  #   avg = calculate_average(score)
  #   std_dev = calculate_std_dev(score.avg)
 #  significance = (float(score) - avg) / std_dev)