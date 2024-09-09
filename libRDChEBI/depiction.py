from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem


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
    useMolBlockWedging=True,
    atomLabelDeuteriumTritium=True
):
    mol = parse_molblock(molfile)
    if not mol:
        return None

    sgs_single_atom = []
    for sg in Chem.GetMolSubstanceGroups(mol):
        sg_props = sg.GetPropsAsDict()
        if sg_props["TYPE"] != "SUP":
            continue
        sg_atoms = list(sg.GetAtoms())
        if len(sg.GetAtoms()) == 1:
            sgs_single_atom.append([sg_atoms, sg_props["LABEL"]])

    for at in mol.GetAtoms():
        dlabel = at.GetSymbol()
        # ChEBI doesn't like to show '#'
        # nor superindices in numbered R groups
        if at.GetAtomicNum() == 0 and len(dlabel) > 1 and dlabel[0] == "R":
            if dlabel[1] == "#":
                at.SetProp("_displayLabel", "R")
            else:
                at.SetProp("_displayLabel", f"R{dlabel[1:]}")
            # add sgroup label if the R group is the only
            # member of a SUP SGROUP
            for sg in sgs_single_atom:
                if at.GetIdx() in sg[0]:
                    at.SetProp("_displayLabel", sg[1])

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
    draw_options.useMolBlockWedging = useMolBlockWedging
    draw_options.atomLabelDeuteriumTritium = atomLabelDeuteriumTritium
    draw.DrawMolecule(mol)
    draw.FinishDrawing()
    svg = draw.GetDrawingText()
    return svg
