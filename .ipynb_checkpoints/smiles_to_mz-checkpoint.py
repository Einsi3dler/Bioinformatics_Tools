import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# 1) Read Excel file
df = pd.read_excel("hmdb.xlsx")

# Convert SMILES column to string to avoid float/NaN issues
df["smiles"] = df["smiles"].astype(str)

def compute_mz_values(smiles):
    """
    Given a SMILES string, return (mz_positive, mz_negative).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    
    mw = Descriptors.ExactMolWt(mol)
    PROTON_MASS = 1.007276
    
    mz_pos = mw + PROTON_MASS  # (M + H)+
    mz_neg = mw - PROTON_MASS  # (M - H)-
    
    return mz_pos, mz_neg

positive_list = []
negative_list = []

for smile in df["smiles"]:
    mz_pos, mz_neg = compute_mz_values(smile)
    positive_list.append(mz_pos)
    negative_list.append(mz_neg)

df["MZ VALUE (Positive)"] = positive_list
df["MZ VALUE (Negative)"] = negative_list

df.to_excel("compounds_with_mz.xlsx", index=False)

print("Done! Updated file: compounds_with_mz.xlsx")
