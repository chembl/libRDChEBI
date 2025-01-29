from typing import Optional, List, Tuple
from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem
from indigo.renderer import IndigoRenderer
from indigo import Indigo

indigo = Indigo()
renderer = IndigoRenderer(indigo)


def depict(
    molfile: str,
    width: int = 300,
    height: int = 300,
    output_format: str = "svg",
    baseFontSize: float = -1,
    fixedFontSize: float = -1,
    minFontSize: float = -1,
    maxFontSize: float = -1,
    useCDKAtomPalette: bool = True,
    explicitMethyl: bool = True,
    scaleBondWidth: bool = False,
    addStereoAnnotation: bool = True,
    useMolBlockWedging: bool = True,
    atomLabelDeuteriumTritium: bool = True,
) -> Optional[bytes | str]:
    """Generate an SVG or PNG depiction of a molecule from a molfile.

    Args:
        molfile: A string containing the molecule data in molfile format
        width: Width of the output image in pixels
        height: Height of the output image in pixels
        output_format: Output format, either 'svg' or 'png'
        baseFontSize: Base font size for atom labels (-1 for auto)
        fixedFontSize: Fixed font size for all labels (-1 for variable)
        minFontSize: Minimum font size for atom labels (-1 for no limit)
        maxFontSize: Maximum font size for atom labels (-1 for no limit)
        useCDKAtomPalette: Use CDK atom colors if True
        explicitMethyl: Show explicit methyl groups if True
        scaleBondWidth: Scale bond widths with drawing size if True
        addStereoAnnotation: Add stereochemistry annotations if True
        useMolBlockWedging: Use molblock wedging info for stereo bonds
        atomLabelDeuteriumTritium: Show D and T labels for deuterium/tritium

    Returns:
        An SVG string or PNG bytes representation of the molecule, or None if parsing fails
    """
    if output_format.lower() not in ["svg", "png"]:
        raise ValueError("Output format must be either 'svg' or 'png'")

    mol = parse_molblock(molfile)
    if not mol:
        return None

    sgs_single_atom: List[Tuple[List[int], str]] = []
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

    drawer_class = rdMolDraw2D.MolDraw2DSVG if output_format.lower() == "svg" else rdMolDraw2D.MolDraw2DCairo
    draw = drawer_class(width, height)
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
    return draw.GetDrawingText()


def depict_indigo(molfile, height=300, width=300, output_format="png"):
    if output_format.lower() not in ["svg", "png"]:
        raise ValueError("Output format must be either 'svg' or 'png'")

    indigo.setOption("render-image-width", width)
    indigo.setOption("render-image-height", height)
    indigo.setOption("render-implicit-hydrogens-visible", True)
    indigo.setOption("render-coloring", True)
    indigo.setOption("render-background-color", "1, 1, 1")
    indigo.setOption("render-output-format", output_format.lower())
    indigo.setOption("render-stereo-style", "none")

    try:
        mol = indigo.loadMolecule(molfile)
        buffer = renderer.renderToBuffer(mol)

        if output_format.lower() == "svg":
            return buffer.decode("utf-8")
        else:  # PNG
            return buffer
    except Exception as e:
        raise Exception(f"Error processing molecule: {str(e)}")
