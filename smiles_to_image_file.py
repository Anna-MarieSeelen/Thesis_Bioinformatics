from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions, MolDrawOptions

def vis_MS2Query_mol(smiles, spectrum_num, motif, substructure_smiles):
    # creating an atom list for visualization
    mol = Chem.MolFromSmiles(smiles)
    substructure_smiles=Chem.MolFromSmiles(smiles)
    mol_block = Chem.MolToMolBlock(mol)
    #atom_list = list(mol.GetSubstructMatch(substructure_smiles))
    atom_list=[0,1,2,3,4,5,6,7,8,9,18]
    bond_list=get_bond_list(atom_list, mol_block)

    opts = MolDrawOptions()
    opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})

    Draw.MolToFile(mol, f"/lustre/BIF/nobackup/seele006/MS2Query_identified_structures/{spectrum_num}_{motif}.png", highlightAtoms=atom_list, highlightBonds=bond_list, options=opts)
    return None

def get_bond_list(atom_list, molblock):
    # creating a bond list for visualization
    mol = Chem.MolFromMolBlock(molblock)

    # creating a bond list for visualization
    bond_list = []
    bonds_in_prec_mol = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
    for bond in bonds_in_prec_mol:
        atom_1,atom_2=bond
        #if atom 1 and atom 2 are also both in the atom_list then they are connected to each other in the molecule!
        if all(x in atom_list for x in [atom_1,atom_2]):
            bond_list.append(mol.GetBondBetweenAtoms(atom_1, atom_2).GetIdx())
    return bond_list


def main():
    print("hiiiiiiiii")
    smiles=argv[1]
    motif=argv[2]
    spectrum_num=argv[3]
    substructure_smiles=argv[4]
    vis_MS2Query_mol(smiles, spectrum_num, motif, substructure_smiles)

if __name__ == "__main__":
    main()