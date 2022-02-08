import os
import sys
import urllib.request
import argparse
import time
from typing import Union, List
from tqdm.auto import tqdm
from joblib import Parallel, delayed


class PDBLoader(object):
    """
    Reading from a protein `PDB IDs` list and download file to specific dir
    """
    def __init__(
        self, 
        filetypes: Union[str, List[str]]='pdb', 
        download_gz: Union[bool, List[bool]]=False, 
        datadir="./data"
    ):
        self.datadir = datadir
        if isinstance(filetypes, str):
            filetypes = [filetypes]
        if isinstance(download_gz, bool):
            download_gz = [download_gz] * len(filetypes)
        url_prefix = []
        file_suffix = []
        for ftype, squeeze in zip(filetypes, download_gz):
            if ftype == 'pdb':
                url_prefix.append("https://files.rcsb.org/download/")
                if squeeze: file_suffix.append(".pdb.gz")
                else: file_suffix.append(".pdb")
            elif ftype == 'pdbx/mmcif':
                url_prefix.append("https://files.rcsb.org/download/")
                if squeeze: file_suffix.append(".cif.gz")
                else: file_suffix.append(".cif")
            elif ftype == 'pdbml/xml':
                url_prefix.append("https://files.rcsb.org/download/")
                file_suffix.append(".xml")
                assert squeeze == True
            elif ftype == 'fasta':
                url_prefix.append("https://www.rcsb.org/fasta/entry/")
                file_suffix.append("")
                assert squeeze == False
            else:
                raise ValueError(f'Invalid type {ftype}.')
        self.filetypes = filetypes
        self.file_suffix = file_suffix
        self.url_prefix = url_prefix
    
    def __call__(self, pdbcode: str):
        assert pdbcode is not None

        for file, object_url in zip(self.file_suffix, self.url_prefix):
            url = object_url + pdbcode + file
            outpath = os.path.join(self.datadir, pdbcode)
            outfnm = os.path.join(outpath, pdbcode + file)
            if not os.path.exists(outpath): os.makedirs(outpath)
            try:
                urllib.request.urlretrieve(url, outfnm)
            except Exception as err:
                print(str(err), file=sys.stderr)
                continue


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_list', type=str)
    parser.add_argument('--filetype', nargs='+')
    parser.add_argument('--load_gz', action='store_true', default=False)
    parser.add_argument('--dir', type=str, default='./data')
    parser.add_argument('--num_jobs', type=int, default=16)
    args = parser.parse_args()

    loader = PDBLoader(filetypes=args.filetype, download_gz=args.load_gz, datadir=args.dir)
    parallel = Parallel(n_jobs=args.num_jobs)
    with open(args.pdb_list, 'r') as pdb_list:
        pdb = pdb_list.read().rstrip('\n').split(', ')
    print('Loading begin...')
    start = time.time()
    parallel(delayed(loader)(protein) for protein in tqdm(pdb))
    cost = start - time.time()
    print(f'Finished in {cost // 60} min {cost % 60:.1f} s!')
