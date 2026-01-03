import random
from rdkit import Chem
from rdkit.Chem import AllChem

# --- 1. HAM MADDELER (BUILDING BLOCKS) ---
# Markette zaten var olduğunu varsaydığımız parçalar
acids = [
    "CC(=O)O", "CCC(=O)O", "c1ccccc1C(=O)O", "CC(C)C(=O)O",
    "O=C(O)C1CCCCC1", "FC(=O)O", "ClCC(=O)O"
]
alcohols = [
    "CO", "CCO", "CCCO", "CC(C)O", "c1ccccc1O",
    "OC1CCCCC1", "FC(F)(F)CO"
]
amines = [
    "CN", "CCN", "CCCN", "c1ccccc1N", "NC1CCCCC1",
    "CC(C)N", "ClCCN"  # DÜZELTİLDİ: "ClcN" -> "ClCCN" (Kloroetilamin)
]
aldehydes = [
    "C=O", "CC=O", "CCC=O", "c1ccccc1C=O", "O=Cc1ccccc1F"
]

# --- 2. SENTEZ REAKSİYONLARI (Bizim Kuralların Tersi) ---
# Ajan parçalamayı öğreniyor, biz burada birleştirmeyi yapıyoruz.

reactions = {
    "Ester": AllChem.ReactionFromSmarts("[C:1](=[O:2])[OH].[O:3][#6:4]>>[C:1](=[O:2])-[O:3]-[#6:4]"),
    # Asit + Alkol -> Ester
    "Amide": AllChem.ReactionFromSmarts("[C:1](=[O:2])[OH].[N:3][#6:4]>>[C:1](=[O:2])-[N:3]-[#6:4]"),
    # Asit + Amin -> Amid
    "Imine": AllChem.ReactionFromSmarts("[C:1]=[O:2].[N:3]>>[C:1]=[N:3]")  # Aldehit + Amin -> İmin
}

generated_molecules = []

print("--- Molekül Üretimi Başlıyor ---\n")

# 1. ESTERLER (Asit + Alkol)
for acid in acids:
    for alc in alcohols:
        mol_acid = Chem.MolFromSmiles(acid)
        mol_alc = Chem.MolFromSmiles(alc)

        # Güvenlik Kontrolü: Eğer SMILES hatalıysa atla
        if mol_acid is None or mol_alc is None:
            continue

        ps = reactions["Ester"].RunReactants((mol_acid, mol_alc))
        if ps:
            prod = ps[0][0]
            try:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                generated_molecules.append({"smiles": smi, "type": "Ester", "parents": f"{acid} + {alc}"})
            except:
                continue

# 2. AMİDLER (Asit + Amin)
for acid in acids:
    for amine in amines:
        mol_acid = Chem.MolFromSmiles(acid)
        mol_amine = Chem.MolFromSmiles(amine)

        # Güvenlik Kontrolü
        if mol_acid is None or mol_amine is None:
            continue

        ps = reactions["Amide"].RunReactants((mol_acid, mol_amine))
        if ps:
            prod = ps[0][0]
            try:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                generated_molecules.append({"smiles": smi, "type": "Amide", "parents": f"{acid} + {amine}"})
            except:
                continue

# 3. İMİNLER (Aldehit + Amin)
for ald in aldehydes:
    for amine in amines:
        mol_ald = Chem.MolFromSmiles(ald)
        mol_amine = Chem.MolFromSmiles(amine)

        # Güvenlik Kontrolü
        if mol_ald is None or mol_amine is None:
            continue

        ps = reactions["Imine"].RunReactants((mol_ald, mol_amine))
        if ps:
            prod = ps[0][0]
            try:
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                generated_molecules.append({"smiles": smi, "type": "Imine", "parents": f"{ald} + {amine}"})
            except:
                continue

# --- SONUÇLARI YAZDIR ---
print(f"Toplam {len(generated_molecules)} adet yeni molekül üretildi!\n")
print("# Bu listeyi mcts_agent.py dosyasındaki TEST_MOLECULES kısmına kopyalayabilirsin.")
print("TEST_MOLECULES = [")

# Rastgele 50 tane seçelim (veya hepsini)
sample_size = min(50, len(generated_molecules))
selected = random.sample(generated_molecules, sample_size)

for i, mol in enumerate(selected):
    print(
        f'    {{"name": "Sentetik Molekül {i + 1}", "smiles": "{mol["smiles"]}", "expected_rule": "{mol["type"]} Reaction"}},')

print("]")