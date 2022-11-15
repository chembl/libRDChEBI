from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem.Draw import rdMolDraw2D


def depict(molfile, width=400, height=400):
    mol = parse_molblock(molfile)
    d = rdMolDraw2D.MolDraw2DSVG(width, height)
    draw_options = d.drawOptions()
    draw_options.minFontSize = 8
    draw_options.maxFontSize = 14
    draw_options.explicitMethyl = True
    draw_options.scaleBondWidth = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    p = d.GetDrawingText()
    return p
