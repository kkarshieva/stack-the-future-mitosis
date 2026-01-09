import os
from dotenv import load_dotenv
from supabase import create_client
from app import save_match

load_dotenv()

SUPABASE_URL = os.environ["SUPABASE_URL"]
SUPABASE_SERVICE_ROLE_KEY = os.environ["SUPABASE_SERVICE_ROLE_KEY"]

sb = create_client(SUPABASE_URL, SUPABASE_SERVICE_ROLE_KEY)

user_id = "4df694c4-898e-4983-8bc3-891b312ef680"
compound_name = "TestCompound"
smiles = "O=C(C)Oc1ccccc1C(=O)O"
properties = {
    "molecular_weight": "101-200",
    "lipophilicity": "1-3",
    "hydrogen_bonding_acceptors": "<5",
    "hydrogen_bonding_donors": ">10",
}

save_match(sb, user_id, compound_name, smiles, properties)
print("Test compound saved successfully!")
