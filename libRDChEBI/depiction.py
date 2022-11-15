from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem.Draw import rdMolDraw2D
import ctfile
import io


def depict(molfile, height=400, width=400):
    mol = parse_molblock(molfile)

    # hack
    f = io.StringIO(molfile)
    ctf = ctfile.load(f)
    r_idxs = []
    for idx, at in enumerate(ctf.atoms):
        if at.atom_symbol in ("R", "R#"):
            r_idxs.append(idx)
    for at in mol.GetAtoms():
        if at.GetIdx() in r_idxs and at.GetSymbol() == "*":
            at.SetProp("atomLabel", "R")
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
