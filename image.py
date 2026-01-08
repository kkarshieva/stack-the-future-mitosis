from rdkit import Chem
from rdkit.Chem import Draw

FILTERED_PATH = "datasets/compounds_filtered.csv"
df = pd.read_csv(FILTERED_PATH)
df = df[["Name", "Smiles"]]  # keeping only these two for images
compounds = df.to_dict(orient="records")
