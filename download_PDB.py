import urllib.request
import urllib.error
def main():
    base_url = "http://files.rcsb.org/view/"
    pdb_id = "4H6Z"
    target_url = base_url+pdb_id.strip()+".pdb"
    pdb_file = "./"+pdb_id.strip()+".pdb"
    try:
        urllib.request.urlretrieve(target_url,pdb_file)
        print(pdb_file+"has been download")
    except(urllib.error.HTTPError or TimeoutError or urllib.error.URLError):
        print("ERROR: Network connection problem when "
              "trying to download " + pdb_file + ". \nReturning to program.." )
main()