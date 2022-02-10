import os
import sys
import socket
import re
import shutil
import urllib.request
import argparse
import time
import itertools
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
                file_suffix.append(".fasta")
                assert squeeze == False
            else:
                raise ValueError(f'Invalid type {ftype}.')
        self.filetypes = filetypes
        self.file_suffix = file_suffix
        self.url_prefix = url_prefix
    
    def __call__(self, pdbcode: str):
        assert pdbcode is not None

        missing = []

        for file, object_url in zip(self.file_suffix, self.url_prefix):
            url = object_url + pdbcode + file
            outpath = os.path.join(self.datadir, pdbcode)
            outfnm = os.path.join(outpath, pdbcode + file)
            if not os.path.exists(outpath): os.makedirs(outpath)
            try:
                urllib.request.urlretrieve(url, outfnm)
            except socket.timeout:
                count = 1
                while count <= 3:
                    try:
                        urllib.request.urlretrieve(url, outfnm)
                        break
                    except socket.timeout:
                        print(f'Reloading for {count} times on {url}', file=sys.stdout)
                        count += 1
                if count > 3:
                    missing.append(pdbcode)
                    print(f'Loading time out for {pdbcode}!')
            except Exception as err:
                print(str(err) + f' when downloading {pdbcode}', file=sys.stderr)
                missing.append(pdbcode)
                continue
        
        return missing


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_list', type=str)
    parser.add_argument('--split', type=str, default=None)
    parser.add_argument('--re', type=str, default=None)
    parser.add_argument('--filetype', nargs='+')
    parser.add_argument('--load_gz', action='store_true', default=False)
    parser.add_argument('--dir', type=str, default='./data')
    parser.add_argument('--num_jobs', type=int, default=16)
    parser.add_argument('--time_out', type=int, default=15)
    args = parser.parse_args()

    loader = PDBLoader(filetypes=args.filetype, download_gz=args.load_gz, datadir=args.dir)
    parallel = Parallel(n_jobs=args.num_jobs)
    socket.setdefaulttimeout(args.time_out)
    with open(args.pdb_list, 'r') as pdb_list:
        pdb = pdb_list.read().rstrip('\n').split(args.split)
    if args.re is not None:
        pdb = [re.findall(args.re, s) for s in pdb]
    print('Loading begin...')
    start = time.time()
    missing_p = parallel(delayed(loader)(protein) for protein in tqdm(pdb))
    cost = start - time.time()
    print(f'Finished in {cost // 60} min {cost % 60:.1f} s!')

    missing = list(itertools.chain(*missing_p))
    if len(missing) > 0:
        print(f'Log missing downloading target... missing target: {len(missing)}!')
        with open(os.path.join(args.dir, 'missing.log')) as log:
            log.write(' '.join(missing))
        while True:
            flag = input('Deleting empty file[Y/n]?\n')
            if flag in ['Y', 'n']: break
            else: 
                print('Invalid input!\b')
                time.sleep(0.5)
        if flag == 'Y':
            for pdbcode in missing:
                shutil.rmtree(os.path.join(args.dir, pdbcode))
