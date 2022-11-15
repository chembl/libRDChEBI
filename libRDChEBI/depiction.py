from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem.Draw import rdMolDraw2D
import ctfile
import io


def depict(molfile, height=400, width=400):
    mol = parse_molblock(molfile)

    # hack to get back R's without 'M RGP' entry in molfile
    f = io.StringIO(molfile)
    ctf = ctfile.load(f)
    for at, rd_at in zip(ctf.atoms, mol.GetAtoms()):
        if at.atom_symbol[0] == "R" and rd_at.GetSymbol() == "*":
            rd_at.SetProp("atomLabel", "R")
    # hack

    d = rdMolDraw2D.MolDraw2DSVG(height, width)
    draw_options = d.drawOptions()
    draw_options.minFontSize = 8
    draw_options.maxFontSize = 14
    draw_options.explicitMethyl = True
    draw_options.scaleBondWidth = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    p = d.GetDrawingText()
    return p
