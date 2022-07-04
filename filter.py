"""
PDB file filter
"""

import argparse
import os
import sys
from tqdm.auto import tqdm
from Bio.PDB import MMCIFParser, PDBParser, Selection
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Polypeptide import three_to_one, is_aa


NON_STANDARD_SUBSTITUTIONS = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}


def augmented_is_aa(three):
    if three in NON_STANDARD_SUBSTITUTIONS:
        three = NON_STANDARD_SUBSTITUTIONS[three]
    return is_aa(three, standard=True)


def augmented_three_to_one(three):
    if three in NON_STANDARD_SUBSTITUTIONS:
        three = NON_STANDARD_SUBSTITUTIONS[three]
    return three_to_one(three)


def augmented_get_sequence(chain):
    sequence_from_pdb = list()
    for res in chain:
        resname = res.get_resname()
        if not augmented_is_aa(resname): continue
        sequence_from_pdb.append(augmented_three_to_one(resname))
    
    return ''.join(sequence_from_pdb)


def chain_split(root, pdb_id):
    parser = PDBParser()
    protein_path = os.path.join(root, pdb_id + '.pdb')
    io = PDBIO()
    if not os.path.exists(protein_path):
        return None
    try:
        structure = parser.get_structure(pdb_id, protein_path)
    except Exception:
        print(f'Error in reading {pdb_id}!', file=sys.stderr)
        return None
    chains = Selection.unfold_entities(structure, 'C')
    chains_dict = {chain.get_id(): chain for chain in chains}
    seq_dict = {chain.get_id(): augmented_get_sequence(chain) for chain in chains}
    for chain_id, select_chain in chains_dict.items():
        save_path = os.path.join(root, pdb_id+'_'+chain_id+'.pdb')
        io.set_structure(select_chain)
        try:
            io.save(save_path)
        except:
            print(f'save error in {pdb_id}')
    for chain_id, select_seq in seq_dict.items():
        save_path = os.path.join(root, pdb_id+'_'+chain_id+'.fasta')
        fasta = f'>{pdb_id}_{chain_id}\n' + select_seq + '\n'
        with open(save_path, 'w') as handle:
            handle.write(fasta)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset', type=str)
    # parser.add_argument('--inplace', action='store_true', default=False)
    args = parser.parse_args()

    root = args.dataset
    for file in tqdm(os.listdir(root)):
        if len(file) == 4:
            chain_split(os.path.join(root, file), file)
