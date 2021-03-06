"""
Walk through two data sets and report outputs
"""
import re
import os
import subprocess
import errno
import logging
import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm
from typing import List

_METHOD_DICT = {}

def register_method(name):
    def decorator(cls):
        _METHOD_DICT[name] = cls
        return cls
    return decorator


def get_methods(backends, namelist):
    return {name: _METHOD_DICT[name](backend) for backend, name in zip(backends, namelist)}


def setup_backend(path=None):
    if path is None:
        raise ValueError(
            'Please specify a path!',
        )
    if os.path.isfile(path):
        return path
    else:
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT),
            path
        )


def setup_dataset(path=None):
    if path is None:
        raise ValueError(
            'Please specify a path!',
        )
    if os.path.exists(path):
        dataset = dict()
        for root, dirs, files in os.walk(path):
            if len(files) == 0: continue
            for file in files:
                if file[-4:] == '.pdb':
                    dataset[file[:-4]] = os.path.join(root, file)
        assert len(dataset) > 0, 'Data path contains NO PDB file!'
        return dataset
    else:
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT),
            path
        )
    

class BenchmarkScore(object):
    def __init__(self, path):
        self.path = setup_backend(path)
        self.output = dict()
    
    def __call__(self, pdb1, pdb2):
        """
        Will use PDB2 as reference if needed
        """
        if os.path.isfile(pdb1) and os.path.isfile(pdb2):
            out = subprocess.check_output([self.path, pdb1, pdb2])
            self.lines = str(out).split('\\n')
        else:
            raise Exception(f'Check if {pdb1} and {pdb2} exits!')


# =============================================================================
# Benchmark methods should be added here
# =============================================================================
@register_method('tmalign')
class TMalign(BenchmarkScore):
    def __init__(self, path):
        super().__init__(path)
        self.tm_score = list()
    
    def __call__(self, pdb1, pdb2):
        super().__call__(pdb1, pdb2)
        for line in self.lines:
            data = re.sub(r'\s\s+', ' ', line).split(' ')
            if data[0] == 'TM-score=':
                self.tm_score.append(float(data[1]))
        
        return self.tm_score[-1], self.lines


@register_method('tmscore')
class TMscore(BenchmarkScore):
    def __init__(self, path):
        super().__init__(path)
    
    def __call__(self, pdb1, pdb2):
        super().__call__(pdb1, pdb2)
        for line in self.lines:
            data = re.sub(r"\s\s+", " ", line).split(' ')
            if data[0] == "TM-score" and data[1] == "=":
                self.tm_score = float(data[2])
        
        return self.tm_score, self.lines
# =============================================================================
# Benchmark methods should be added here
# =============================================================================


class Benchmark(object):
    def __init__(self, obj_path, ref_path, log, backends: List[str] = None, methods: List[str] = None):
        """
        Args:
            obj_path: object path that contains all object pdb
            ref_path: reference path that contains all reference pdb
            All matching files should be like `xxx.pdb` with the same file name
        """
        self.obj_dataset = setup_dataset(obj_path)
        self.ref_dataset = setup_dataset(ref_path)
        # find matching keys
        self.matching_datapoints = list(
            self.obj_dataset.keys() & self.ref_dataset.keys())
        assert len(methods) > 0, 'Please specify at least one method'
        self.methods = get_methods(backends, methods)
        if not os.path.exists(log): os.mkdir(log)
        self.log_path = log
        logging.basicConfig(
            filename=os.path.join(log, 'output.log'), encoding='utf-8')

    def __call__(self, num_jobs=128):
        score_frame = dict()
        parallel = Parallel(n_jobs=num_jobs)
        for method in self.methods:
            print(f'Benchmark {method} begin...')
            logging.info(f'Benchmark for {method}\n# '+ '='*20 + '\n')
            scores, outputs = zip(*parallel(
                delayed(self.methods[method])(
                self.obj_dataset[point], self.ref_dataset[point])
                for point in tqdm(self.matching_datapoints)
            ))
            logging.info(outputs)
            score_frame[method] = scores
        
        return pd.DataFrame(
            score_frame, index=self.matching_datapoints
        )
                

if __name__ == '__main__':
    import argparse
    import pickle

    parser = argparse.ArgumentParser()
    parser.add_argument('object_dataset', type=str)
    parser.add_argument('reference_dataset', type=str)
    parser.add_argument('--backend', nargs='+')
    parser.add_argument('--methods', nargs='+')
    parser.add_argument('--log_dir', type=str, default='./log')
    parser.add_argument('--num_workers', type=int, default=128)
    args = parser.parse_args()

    benchmark = Benchmark(
        args.object_dataset, args.reference_dataset, log=args.log_dir, 
        backends=args.backend, methods=args.methods, )
    score_pd = benchmark(args.num_workers)
    with open(os.path.join(args.log_dir, 'benchmark.pkl'), 'wb') as handle:
        pickle.dump(score_pd, handle, protocol=pickle.HIGHEST_PROTOCOL)
