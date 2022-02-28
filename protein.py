"""
Augmented protein metadata

Fasta file will converted into Fasta class in a easy to use dict form

Protein level:
pdb id
    protein pdb id
chain id
    protein chain id

Sequence level:
amino acid
    amino acid type
meta data
    This meta data contains 
"""


import warnings
import re
from easydict import EasyDict
from Bio import BiopythonWarning
from Bio import Seq, SeqRecord
from typing import List, Optional


class FastaData(object):
    def __init__(
        self,
        pdb_id: Optional[str]=None,
        chain_id: Optional[List[str]]=None,
        sequence: Optional[Seq.Seq]=None,
        meta: Optional[dict]=None,
    ):
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.sequence = sequence
        if meta is not None:
            self.set_meta(meta)
    
    @property
    def name(self):
        return f'{self.pdb_id}_{self.chain_id[0]}'
    
    def set_meta(self, meta_info):
        assert isinstance(meta_info, dict)
        self.meta = EasyDict(meta_info)

    
# =============================================================================
# Reading methods should be added here
# =============================================================================
    @classmethod
    def from_pdb_fasta(cls, record):
        """
        Read form a fasta SeqIO record <class 'Bio.SeqRecord.SeqRecord'>
        """
        warnings.simplefilter('ignore', BiopythonWarning)
        assert isinstance(record, SeqRecord.SeqRecord)
        pdb_id = re.match(r'[1-9][A-Za-z0-9]{3}[^A-Za-z0-9]', record.id)
        name, chain_id = None, None
        if pdb_id is not None:
            pdb_id = pdb_id.group(0)[:-1].upper()
            name = re.search(pdb_id + r'_[A-Za-z0-9]', record.name)
        if name is not None:
            name = name.group(0)
        if record.description is not None:
            chain_id = re.search(r'[Cc]hains? [A-Z]( ,?[A-Z])*', record.description)
        if chain_id is not None:
            chain_id = re.sub(r'[Cc]hains?', '', chain_id.group(0))
            chain_id = list(chain_id.replace(' ', '').replace(',', ''))
        assert record.seq is not None
        meta = {'name': name, 'description': record.description}

        return cls(pdb_id, chain_id, record.seq, meta)
# =============================================================================
# Reading methods should be added here
# =============================================================================

