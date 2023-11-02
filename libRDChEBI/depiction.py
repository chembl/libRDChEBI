from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem.Draw import rdMolDraw2D


def depict(
    molfile,
    width=300,
    height=300,
    baseFontSize=-1,
    fixedFontSize=-1,
    minFontSize=-1,
    maxFontSize=-1,
    useCDKAtomPalette=True,
    explicitMethyl=True,
    scaleBondWidth=False,
    addStereoAnnotation=True,
):
    mol = parse_molblock(molfile)
    if not mol:
        return None

    # ChEBI doesn't like to show '#'
    # nor superindices in numbered R groups
    for at in mol.GetAtoms():
        dlabel = at.GetSymbol()
        if len(dlabel) > 1 and dlabel[0] == "R":
            if dlabel[1] == "#":
                at.SetProp("_displayLabel", "R")
            else:
                at.SetProp("_displayLabel", f"R{dlabel[1:]}")

    draw = rdMolDraw2D.MolDraw2DSVG(width, height)
    draw_options = draw.drawOptions()
    draw_options.baseFontSize = baseFontSize
    draw_options.fixedFontSize = fixedFontSize
    draw_options.useCDKAtomPalette = useCDKAtomPalette
    draw_options.minFontSize = minFontSize
    draw_options.maxFontSize = maxFontSize
    draw_options.explicitMethyl = explicitMethyl
    draw_options.scaleBondWidth = scaleBondWidth
    draw_options.addStereoAnnotation = addStereoAnnotation
    draw.DrawMolecule(mol)
    draw.FinishDrawing()
    svg = draw.GetDrawingText()
    return svg
