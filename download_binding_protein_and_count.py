import urllib.request
import re
import os
import operator
pdb_list = {"1JFF", "2HXF", "2HXH", "2P4N", "2WBE", "3DCO", "3EDL", "3IZ0", "3J6G", "3J6P" ,"4AQV", "4AQW", "4CK5",
          "4CK6", "4CK7", "4UXO", "4UXP", "4UXR", "4UXS", "4UXT", "4UXY", "4UY0", "5HNW", "5HNX", "5HNY", "5HNZ",
          "5M50", "5M54", "5M5C", "5M5I", "5M5L", "5M5M", "5M5N", "5M5O", "5ND2", "5ND3", "5ND4", "5ND7", "5SYE", "5SYF", "5W3J"}
thr = 4.
def download_pdb(pdb_id):
    base_url = "https://files.rcsb.org/view/"
    file_name = pdb_id + ".pdb"
    target_url = base_url + file_name
    directory_path = "/Users/leilou/Desktop/"
    pdb_file = directory_path + file_name
    if (file_name.lower() in os.listdir(directory_path)) or (
                file_name.upper() in os.listdir(directory_path)):
        print("\n" + "-----------------\n" +
              file_name + " is already in " + directory_path)  # +
    else:
        try:
            urllib.request.urlretrieve(target_url, pdb_file)
            print(file_name + " downloaded>>>>")
            #input("Pause01! Press enter to start again.")
        except(urllib.error.HTTPError or TimeoutError or urllib.error.URLError):
            print("ERROR: Network connection problem when "
                  "trying to download " + pdb_file + " from " + target_url + ". \nReturning to program...")

def parse_pdb(pdb_id):
    file_name = pdb_id + ".pdb"
    directory_path = "/Users/leilou/Desktop/"
    pdb_file = directory_path + file_name
    atom = []
    het  = []
    res = []
    ll = []
    na = 0 #no. of atom
    nh = 0 #no. of het atom
    nr = 0 #no. of contact residues
    caa = {} # contact amino acid
    aap = {"ALA": 0,"ARG": 0,"CYS": 0,"ASP": 0,"ASN": 0,"GLU": 0,"GLN": 0,"MET": 0,"LYS": 0,"HIS": 0,"SER": 0,"THR": 0,"PHE": 0,"TYR": 0,"VAL": 0,"GLY": 0,"LEU": 0,"ILE": 0,"TRP": 0,"PRO": 0,"ALL": 0,} # amino acid preferences
    with open(pdb_file) as f:
        for line in f:
            if re.match('.*RESIDUE.*',line):
                 res.append(" ")
                 res[nr] = [" "]*6
                 res[nr][0] = line[17:20].strip()+'-'+line[22:26].strip()+'-'+line[21]
                 res[nr][1] = line[30:38].strip()  # x
                 res[nr][2] = line[38:46].strip()  # y
                 res[nr][3] = line[46:54].strip()  # z
                 res[nr][4] = line[17:20].strip()  # res name
                 res[nr][5] = line[76:78].strip()  # element
                 print(line)
            if re.match('^ATOM.*', line):
                atom.append(" ")
                atom[na] = [" "]*6
                atom[na][0] = line[17:20].strip()+'-'+line[22:26].strip()+'-'+line[21]
                atom[na][1] = line[30:38].strip()  # x
                atom[na][2] = line[38:46].strip()  # y
                atom[na][3] = line[46:54].strip()  # z
                atom[na][4] = line[17:20].strip()  # res name
                atom[na][5] = line[76:78].strip()  # element

                na += 1
            elif re.match('^HETATM.*TA1', line):
                het.append(" ")
                het[nh] = [" "] * 4
                het[nh][0] = line[17:20].strip()+'-'+line[22:26].strip()+'-'+line[21]+'-'+line[12:16].strip() #res name+res no.+chain ID
                het[nh][1] = line[30:38].strip() # x
                het[nh][2] = line[38:46].strip() # y
                het[nh][3] = line[46:54].strip() # z
                nh += 1
            elif re.match('^(HEADER|TITLE)',line):
                print(line)
                #input("pause01: press Enter to continue")
    #print("number of atoms:"+str(na)+" number of het atoms:"+str(nh))
    for i in range(na):
        if atom[i][0] in caa:
            continue #
        for j in range(nh):
            dist = ((float(atom[i][1])-float(het[j][1]))**2+(float(atom[i][2])-float(het[j][2]))**2+(float(atom[i][3])-float(het[j][3]))**2)**.5
            if dist < thr:
                #print(dist,atom[i][0],het[j][0])
                caa[atom[i][0]] = 1
                nr += 1
                aap[atom[i][4]] += 1
                aap['ALL']+=1
                break
    print("no. of contact residues:"+str(nr))
    for i in aap.keys():
        print (i, aap[i], str(aap[i] / aap['ALL'] * 100.))
    sorted_x = sorted(aap.items(), key=lambda x: x[1], reverse=True)
    print("The most enriched 3 AA are: ")
    print(sorted_x[1:4])

def main():
    for pdb_id in pdb_list:
        if re.match('^\w{4}$',pdb_id):
            download_pdb(pdb_id)
            parse_pdb(pdb_id)
main()
