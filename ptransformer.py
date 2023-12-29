'''
Transform anything in protein to make things simple
'''

import os
import warnings

import numpy as np
from tmtools import tm_align
from tmtools.io import get_residue_data
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Polypeptide import three_to_one

STD_AA = ['GLY', 'ALA', 'SER', 'THR', 'ASP', 'GLU', 'VAL', 'ILE', 
          'TYR', 'TRP', 'LEU', 'PHE', 'ASN', 'GLN', 'PRO', 'MET', 
          'CYS', 'LYS', 'ARG', 'HIS']

def isRotationMatrix(M: np.ndarray):
    tag = False
    I = np.identity(M.shape[0])
    if np.all((np.matmul(M, M.T)) == I) and (np.linalg.det(M)==1): tag = True
    return tag   

class TransFormer(object):
    """
        Tools to rename PDB chains, trim PDB structure and align sequence and structure
    """
    def __init__(self, pdb_id: str, pdb_path: str, index_gap: int=32) -> None:
        """
        Args:
            pdb_id (str): pdb id of the pdb file you want to deal with
            pdb_path (str): path of the pdb file
            index_gap (int, optional): index gap that will simulate a chain split
        """
        self.pdb_id = pdb_id
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure(self.pdb_id, pdb_path)
        self.index_gap = index_gap
        
    def get_chain_seq(self, chain_id: str) -> str:
        """
            get sequence with NO NON-STANDARD amino acids
        Args:
            chain_id (str): the chain id you want to get sequence
        """
        res_list = list()
        for res in self.structure[0][chain_id]:
            if res.get_resname() in STD_AA:
                res_list.append(three_to_one(res.get_resname()))
                
        return ''.join(res_list)
        
    def rename_chain(self, rename_map: dict, update: bool=False) -> Structure:
        """
            rename PDB chain using the mapping of rename_map and return a standard structure
        Args:
            rename_map (dict): Name mapping, if the chain id in PDB file is not provided,
            the name will keep the same. If none of the chain id can match the key in name
            mapping, there will be a warning. You can provide a mapping that map different
            chain to same id, this will ends up with a chain merge with a default 32 index gap.
            update (bool, optional): Whether to update self.structure. Defaults to False.

        Returns:
            Structure: A structure of Biopython. Rename the chain id with also the index and
            non-standard amino acids
        """
        builder = StructureBuilder()
        builder.init_structure(self.pdb_id)
        builder.init_model(0)
        index_dict = {key: 0 for key in rename_map}
        def _build_chain(chain, chain_id):
            if index_dict[chain_id] != 0:
                index_dict[chain_id] += self.index_gap
            builder.init_chain(chain_id)
            builder.init_seg('    ')
            for res in chain:
                if res.get_resname() in STD_AA:
                    index_dict[chain_id] += 1
                    builder.init_residue(
                        resname=res.get_resname(),
                        field=' ',
                        resseq=index_dict[chain_id],
                        icode=' ',
                    )
                    for atom in res:
                        builder.init_atom(
                            atom.get_name(),
                            atom.get_coord(),
                            atom.get_bfactor(),
                            atom.get_occupancy(),
                            atom.get_altloc(),
                            atom.get_fullname()
                        )
        renamed_chains = list()
        for chain in self.structure[0]:
            if chain.id.strip() in rename_map:
                _build_chain(chain, rename_map[chain.id.strip()])
                renamed_chains.append(chain.id.strip())
            else:
                _build_chain(chain, chain.id.strip())
        
        if len(renamed_chains) < len(rename_map.keys()):
            ignore_chains = set(rename_map.keys()) - set(renamed_chains)
            warnings.warn(f'chain {ignore_chains} are ignored.')
            
        if update:
            self.structure = builder.get_structure()
            
        return builder.get_structure()
    
    def trim_chain(self, trim_dict: dict, update: bool=False) -> Structure:
        """
            trim the structure and reindex it and return a standard structure
        Args:
            trim_dict (dict): Trim mapping like {chain_id: trimed_sequence}, the chain id must
            be in the structure. Any chain not appear in trim_dict will be discard.
            update (bool, optional): Whether to update self.structure.. Defaults to False.

        Returns:
            Structure: A structure of Biopython. Rename the chain id with also the index and
            non-standard amino acids
        """
        builder = StructureBuilder()
        builder.init_structure(self.pdb_id)
        builder.init_model(0)
        def _build_trimed_chain(chain, chain_id, seq, gap):
            builder.init_chain(chain_id)
            builder.init_seg('    ')
            res_seq_before = 0
            res_seq_after = 0
            for res in chain:
                if res_seq_before < gap:
                    res_seq_before += 1
                    continue
                if res.get_resname() in STD_AA:
                    assert three_to_one(res.get_resname()) == seq[res_seq_after]
                    res_seq_after += 1
                    builder.init_residue(
                        resname=res.get_resname(),
                        field=' ',
                        resseq=res_seq_after,
                        icode=' ',
                    )
                    for atom in res:
                        builder.init_atom(
                            atom.get_name(),
                            atom.get_coord(),
                            atom.get_bfactor(),
                            atom.get_occupancy(),
                            atom.get_altloc(),
                            atom.get_fullname()
                        )
        for chain_id in trim_dict:
            structure_seq = self.get_chain_seq(chain_id)
            assert trim_dict[chain_id] in structure_seq, f'chain {chain_id} seq not in structure'
            gap = structure_seq.find(trim_dict[chain_id])
            _build_trimed_chain(self.structure[0][chain_id], chain_id, trim_dict[chain_id], gap)

        if update:
            self.structure = builder.get_structure()
            
        return builder.get_structure()
    
    def rotate_structure(self, rotation_matrix: np.ndarray, bias: np.ndarray, update: bool=False) -> Structure:
        """
            Transform the 3D structure based on rotation matrix and bias
        Args:
            rotation_matrix (np.ndarray): Rotation matrix that will apply to the structure
            bias (np.ndarray): Bias that will apply to the structure
            update (bool, optional): Whether to update self.structure. Defaults to False.

        Raises:
            AssertionError: We will always check if the rotation matrix is valid

        Returns:
            Structure: A structure of Biopython.
        """
        assert np.shape(rotation_matrix) == (3, 3), f'rotation matrix {rotation_matrix} should have shape (3, 3) but {rotation_matrix.shape}'
        if not isRotationMatrix(rotation_matrix):
            raise AssertionError(f'rotation matrix {rotation_matrix} is not a rotation matrix')
        assert np.shape(bias) == (3, ) or np.shape(bias) == (3, 1), f'bias matrix {bias} must have shape (3, ) or (3, 1) but {bias.shape}'
        if np.shape(bias) == (3, 1):
            bias = bias.squeeze()
        
        builder = StructureBuilder()
        builder.init_structure(self.pdb_id)
        builder.init_model(0)
        def _rotate_chain(chain, chain_id):
            builder.init_chain(chain_id)
            builder.init_seg('    ')
            res_seq = 0
            for res in chain:
                if res.get_resname() in STD_AA:
                    res_seq += 1
                    builder.init_residue(
                        resname=res.get_resname(),
                        field=' ',
                        resseq=res_seq,
                        icode=' ',
                    )
                    for atom in res:
                        builder.init_atom(
                            atom.get_name(),
                            np.matmul(atom.get_coord(), rotation_matrix) + bias,
                            atom.get_bfactor(),
                            atom.get_occupancy(),
                            atom.get_altloc(),
                            atom.get_fullname()
                        )
        for chain in self.structure[0]:
            _rotate_chain(chain, chain.id.strip())
            
        if update:
            self.structure = builder.get_structure()
            
        return builder.get_structure()
    
    def align_3D_to(self, chain_id: str, object_3D: Structure, object_chain_id: str, update: bool=False) -> Structure:
        """
            align structure to object 3D coordinates and return the rotated standard structure
        Args:
            chain_id (str): The chain you want to choose to be aligned.
            object_3D (Structure): Bio.PDB.Structure.Structure instance for alignment.
            object_chain_id (str): The object chain id you want to align to. Must be in object_3D.
            update (bool, optional): Whether to update self.structure. Defaults to False.

        Returns:
            Structure: A structure of Biopython.
        """
        ori_coords, ori_seq = get_residue_data(self.structure[0][chain_id])
        obj_coords, obj_seq = get_residue_data(object_3D[0][object_chain_id])
        res = tm_align(ori_coords, obj_coords, ori_seq, obj_seq)
        if update:
            self.rotate_structure(res.u, res.t, update=True)
            return self.structure
        
        return self.rotate_structure(res.u, res.t)
    
    def save_pdb(self, saved_path: str) -> None:
        """
        Args:
            saved_path (str): The path to save pdb file.
        """
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(saved_path)
        
        
class ABTransFormer(TransFormer):
    """
        Transform antibodies.
    """
    def __init__(self, pdb_id: str, pdb_path: str, heavy_chain_id: str, light_chain_id: str, index_gap: int=32) -> None:
        """
        Args:
            pdb_id (str): pdb id of the pdb file you want to deal with
            pdb_path (str): path of the pdb file
            heavy_chain_id (str): heavy chain id
            light_chain_id (str): light chain id
            index_gap (int, optional): index gap that will simulate a chain split
        """
        super.__init__(pdb_id, pdb_path, index_gap)
        self.rename_chain(rename_map={heavy_chain_id: 'H', light_chain_id: 'L'}, update=True)
    
    def number_antibody(self, numbering_scheme: str='IMGT') -> Structure:
        raise NotImplementedError
    
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument()