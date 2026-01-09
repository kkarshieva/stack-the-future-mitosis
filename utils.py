import numpy as np
import warnings
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore", category=DeprecationWarning)

#turn smiles into fingeprint vector 
def smiles_to_vector(smiles: str, n_bits: int = 2048):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits) #outputs bit vector
        #convert to numpy array
        arr = np.zeros((n_bits,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    except Exception:
        return None

#rocchio feedback algorithm
def rocchio_update(current_vec, liked_vectors, disliked_vectors, alpha=0.9, beta=0.5, gamma=0.3):
    if len(liked_vectors) > 0: #if there are likes compute average
        liked_centroid = np.mean(np.vstack(liked_vectors), axis=0)
    else:
        #if no likes use zeros
        liked_centroid = np.zeros_like(current_vec)

    if len(disliked_vectors) > 0:
        disliked_centroid = np.mean(np.vstack(disliked_vectors), axis=0)
    else:
        disliked_centroid = np.zeros_like(current_vec)
    #apply rocchio formula
    return (alpha * current_vec) + (beta * liked_centroid) - (gamma * disliked_centroid)

# #general explain: Every molecule becomes a fingerprint vector. When you like something, your preference vector moves toward it.”
#When you dislike something, it moves away.Over time, the system learns what ‘type’ of molecules you like.”