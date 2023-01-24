import rdkit

rdkversion = rdkit.__version__.split(".")
if rdkversion < ["2022", "09", "4"]:
    raise ValueError("need an RDKit version >= 2022.09.4")
