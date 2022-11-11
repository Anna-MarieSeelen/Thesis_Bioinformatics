

from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions, MolDrawOptions

def vis_substructure_in_precursor_mol(mol_block):
    opts = MolDrawOptions()
    opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})

    mol = Chem.MolFromMolBlock(mol_block)
    Draw.MolToFile(mol, f"/lustre/BIF/nobackup/seele006/MAGMa_illustrations_of_substructures/{identifier}_{motif}.png", options=opts)
    return None

def main():
    path_to_store_images=argv[1]
    list_of_smiles_to_covert_to_img=argv[2]
    motif=argv[3]


if __name__ == "main":
    main()