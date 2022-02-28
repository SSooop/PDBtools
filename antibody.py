"""
Antibody annotation
"""

from Bio import Seq

from .protein import FastaData


class IgFastaData(FastaData):
    def __init__(
        self,
        pdb_id=None,
        chain_id=None,
        sequence=None,
        annotation=None,
        meta=None,
    ):
        super().__init__(pdb_id, chain_id, sequence, meta)
        if annotation is not None:
            assert isinstance(annotation, Seq.Seq)
            self.annotation = annotation
        
    def annotation(
        self,
        path: str,
        method: str = 'paratome',
    ):
        """
        Annotate antibody sequence using supported method

        Args
            path: path to executable file
            method: annotation mathods including
                'kabat': https://en.wikipedia.org/wiki/Kabat_numbering_scheme
                'chothia': http://www.chemogenomix.com/chothia-antibody-numbering
                'imgt': https://www.imgt.org/IMGTindex/numbering.php
                'paratome': https://openebench.bsc.es/tool/paratome
        
        This method will directly call the executable file using the
        object path
        """

        