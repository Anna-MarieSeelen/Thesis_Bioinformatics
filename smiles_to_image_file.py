from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions, MolDrawOptions

def vis_MS2Query_mol(smiles, spectrum_num, motif):
    opts = MolDrawOptions()
    opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, f"/lustre/BIF/nobackup/seele006/MS2Query_identified_structures/{spectrum_num}_{motif}.png", options=opts)
    return None

def main():
    print("hiiiiiiiii")
    smiles=argv[1]
    motif=argv[2]
    spectrum_num=argv[3]
    vis_MS2Query_mol(smiles, spectrum_num, motif)

if __name__ == "__main__":
    main()