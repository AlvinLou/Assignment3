import os
import re
import random
import math

def parse_pdb(pdb_id):
    file_name = pdb_id
    directory_path = "/Users/leilou/Desktop/"
    pdb_file = directory_path + file_name + ".pdb"
    atom = []
    het  = []
    na = 0 #no. of atom
    nh = 0 #no. of het atom
    nr = 0 #no. of contact residues
    caa = {} # contact amino acid
    #aap = {"ALA": 0,"ARG": 0,"CYS": 0,"ASP": 0,"ASN": 0,"GLU": 0,"GLN": 0,"MET": 0,"LYS": 0,"HIS": 0,"SER": 0,"THR": 0,"PHE": 0,"TYR": 0,"VAL": 0,"GLY": 0,"LEU": 0,"ILE": 0,"TRP": 0,"PRO": 0,"ALL": 0,} # amino acid preferences
    with open(pdb_file) as f:
        for line in f:
            if re.match('^ATOM.*', line):
                atom.append(" ")
                atom[na] = [" "]*7
                atom[na][0] = float(line[30:38].strip()) # x
                atom[na][1] = float(line[38:46].strip()) # y
                atom[na][2] = float(line[46:54].strip()) # z
                atom[na][3] = line[12:16]  #
                atom[na][4] = line[21]  # res name
                atom[na][5] = line[17:20] # element
                atom[na][6] = int(line[22:26])
                na += 1
                pass
            if re.match('^HETATM.*TA1', line):
                #if line[21] == 'B':
                 #   continue
                het.append(" ")
                het[nh] = [" "] * 5
                het[nh][0] = float(line[30:38].strip()) #res name+res no.+chain ID
                het[nh][1] = float(line[38:46].strip())# x
                het[nh][2] = float(line[46:54].strip()) # y
                het[nh][3] = line[12:16].strip()# z
                het[nh][4] = line[21] #atom name
                nh += 1
            elif re.match('^(HEADER|COMPND|TITLE)',line):
                #print(line)
                #input("pause01: press Enter to continue")
                pass
    #print("number of atoms:"+str(na)+" number of het atoms:"+str(nh))
    #for i in range(len(het)):
    #    print(het[i][0],het[i][1],het[i][2],het[i][3],het[i][4])
    return atom,het
def random_n(n):
    a = random.uniform(-1.0*n,n)
    return a
def initialization(c,m):
    #c = [[0.]*3 for i in range(4)]
    ix = random_n(m)  # delta x
    iy = random_n(m)  # delta y
    iz = random_n(m)  # delta z
    for i in range(len(c)):
        c[i][0] += ix
        c[i][1] += iy
        c[i][2] += iz
    return c

def translation(c,m):
    dx = random_n(m)    #delta x
    dy = random_n(m)    #delta y
    dz = random_n(m)    #delta z
    for i in range(len(c)):
        c[i][0] += dx
        c[i][1] += dy
        c[i][2] += dz
    return c

def geometry_centet(c,m):
    v = [0., 0., 0.]
    for i in range(len(c)):
        v[0] += c[i][0] / float(len(c))
        v[1] += c[i][1] / float(len(c))
        v[2] += c[i][2] / float(len(c))
        if v[0] > m:
            v[0] = v[0] - 2 * m
        if v[0] < -m:
            v[0] = 2 * m - v[0]
        if v[1] > m:
            v[1] = v[1] - 2*m
        if v[1] < -m:
            v[1] = 2 * m - v[1]
        if v[2] > m:
            v[2] = v[2] - 2*m
        if v[2] < -m:
            v[2] = 2 * m - v[2]
    return c
def rotation(c,m):
    pi = math.pi
    a = random_n(m) / 180. * pi
    b = random_n(m) / 180. * pi
    g = random_n(m) / 180. * pi
    r = [[0]*3 for i in range(3)]
    r[0][0] = math.cos(a) * math.cos(b)
    r[0][1] = math.cos(a) * math.sin(b) * math.sin(g) - math.sin(a) * math.cos(g)
    r[0][2] = math.cos(a) * math.sin(b) * math.cos(g) + math.sin(a) * math.sin(g)
    r[1][0] = math.sin(a) * math.cos(b)
    r[1][1] = math.sin(a) * math.sin(b) * math.sin(g) + math.cos(a) * math.cos(g)
    r[1][2] = math.sin(a) * math.sin(b) * math.cos(g) - math.cos(a) * math.sin(g)
    r[2][0] = math.sin(b) * -1
    r[2][1] = math.cos(b) * math.sin(g)
    r[2][2] = math.cos(b) * math.cos(g)

    n = [[0] * 3 for _ in range(len(c))]
    o = [[0] * 3 for _ in range(len(c))]
    ori = geometry_centet(c,m)
    for i in range(len(c)):
        n[i][0] = c[i][0] - ori[0][0]
        n[i][1] = c[i][1] - ori[0][1]
        n[i][2] = c[i][2] - ori[0][2]
    for i in range(len(c)):
        for j in range(3):
            o[i][j] = n[i][j]

    for i in range(len(c)):
        n[i][0] = r[0][0] * o[i][0] + r[0][1] * o[i][1] + r[0][2] * o[i][2]
        n[i][1] = r[1][0] * o[i][0] + r[1][1] * o[i][1] + r[1][2] * o[i][2]
        n[i][2] = r[2][0] * o[i][0] + r[2][1] * o[i][1] + r[2][2] * o[i][2]

    for i in range(len(c)):
        c[i][0] = n[i][0] + ori[0][0]
        c[i][1] = n[i][1] + ori[0][1]
        c[i][2] = n[i][2] + ori[0][2]

    return c


def main():
    side_length = 100.
    fout = open("TA1.pdb",'w')
    atype = ['C','N','O','S']
    a,c = parse_pdb("5SYE")
    #print(a)
    c = initialization(c,side_length/2)
    for i in range(20):
        for j in range(len(a)):
            if a[j][4] != a[j-1][4] and j > 0:
                fout.write("TER\n")
            line = "{0:6s}{1:5d} {2:4s} {3:3s} {4:1s}{5:4d}    {6:8.3f}{7:8.3f}{8:8.3f}\n".format('ATOM',j + 1, a[j][3], a[j][5], a[j][4], a[j][6], a[j][0], a[j][1], a[j][2])
            #print(line)
            fout.write(line)
        fout.write("TER\n")
        for j in range(len(c)):
            line = "{0:6s}{1:5d} {2:4s} {3:3s} {4:1s}{5:4d}    {6:8.3f}{7:8.3f}{8:8.3f}\n".format('HETATM', j + 1, c[j][3], 'TA1', c[j][4], 1, c[j][0], c[j][1], c[j][2])
            #print(line)
            fout.write(line)
        fout.write("TER\n")
        fout.write("END\n")
        c = rotation(c,60.)
        c = translation(c,1.)
        #print(c)
    fout.close()
main()
