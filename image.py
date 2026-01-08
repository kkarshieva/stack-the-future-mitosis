from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO


def smiles_to_svg(smiles: str, width: int = 400, height: int = 400) -> bytes:
    """
    makes an SVG image of a molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES")

    Chem.rdCoordGen.AddCoords(mol)
    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    # set drawing options on drawer.getOptions()
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()
