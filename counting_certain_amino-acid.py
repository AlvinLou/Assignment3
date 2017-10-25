import urllib.request
import urllib.error
def parse_pdb(pdb_file):
    with open(pdb_file) as f:
        nof_leu = 0
        for line in f:
            if (line[17:20] == "LEU"):
                nof_leu += 1
        return nof_leu
def main():
    with open("microtubules_list.txt") as f:
        for line in f:
            pdb_list = line.split(", ")
    print(pdb_list)
    nof_leu = [len(pdb_list)]
    n = 0
    for pdb_id in pdb_list:
        base_url = "http://files.rcsb.org/view/"
        target_url = base_url+pdb_id.strip()+".pdb"
        pdb_file = "./"+pdb_id.strip()+".pdb"
        try:
            urllib.request.urlretrieve(target_url,pdb_file)
            print(pdb_file+" has been download")
        except(urllib.error.HTTPError or TimeoutError or urllib.error.URLError):
            print("ERROR: Network connection problem when "
                  "trying to download " + pdb_file + ". \nReturning to program.." )
        nof_leu[n] = parse_pdb(pdb_file)
        print(pdb_file+" has "+str(nof_leu[n])+" LEU")
main()