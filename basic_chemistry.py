import sys

# RDKit kütüphanelerini içe aktarma
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("HATA: RDKit kütüphanesi bulunamadı.")
    print("Lütfen terminalden 'pip install rdkit' komutunu çalıştırın.")
    sys.exit(1)


def main():
    print("--- ADIM 1: Molekül Tanımlama ---")

    # 1. Aspirin'in SMILES kodu
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"

    mol = Chem.MolFromSmiles(aspirin_smiles)

    if mol:
        print(f"Molekül Başarıyla Yüklendi: Aspirin")
        print(f"SMILES Kodu: {aspirin_smiles}")
        print(f"Atom Sayısı: {mol.GetNumAtoms()}")
    else:
        print("Molekül yüklenemedi!")
        return

    print("\n--- ADIM 2: Retrosentez Kuralı (Aksiyon) Tanımlama ---")

    # 2. Reaksiyon Kuralı (SMARTS) - GÜNCELLENDİ
    # Hata Nedeni: [C:4] sadece alifatik (düz) karbonları arıyordu. Aspirin'de halka var (aromatik).
    # Çözüm: [#6:4] kullandık. #6, atom numarası 6 olan (hem aromatik hem alifatik) tüm karbonları kapsar.

    # Kural: [C:1](=[O:2])-[O:3]-[#6:4]  >>  [C:1](=[O:2])O.[#6:4][O:3]
    # Anlamı: Karbonil(1)-Oksijen(3) bağını kopar. Oksijen(3) diğer karbonla(4) kalsın.

    retro_template = "[C:1](=[O:2])-[O:3]-[#6:4]>>[C:1](=[O:2])O.[#6:4][O:3]"

    print(f"Uygulanacak Kural (SMARTS): {retro_template}")
    print("Anlamı: Ester bağını kır -> Asit + Alkol/Fenol elde et.")
    print("(Not: [#6] etiketi sayesinde hem düz zincir hem de halkalı yapılar tanınacak.)")

    try:
        rxn = AllChem.ReactionFromSmarts(retro_template)
    except Exception as e:
        print(f"Kural hatası: {e}")
        return

    print("\n--- ADIM 3: Reaksiyonu Uygulama ---")

    products = rxn.RunReactants((mol,))

    if not products:
        print("Bu kural bu moleküle uygulanamadı! (Eşleşen yapı yok)")
    else:
        print(f"İşlem Başarılı! {len(products)} farklı sonuç bulundu.")

        for i, product_set in enumerate(products):
            print(f"\nSonuç {i + 1} (Olası Parçalanma Yolu):")
            for j, fragment in enumerate(product_set):
                try:
                    # RDKit bazen reaksiyon sonrası molekülü 'kirli' bırakır, temizliyoruz
                    Chem.SanitizeMol(fragment)
                    frag_smiles = Chem.MolToSmiles(fragment)
                    print(f"  - Parça {j + 1}: {frag_smiles}")
                except Exception as e:
                    print(f"  - Parça {j + 1}: {Chem.MolToSmiles(fragment)} (Ham Hali - Temizlenemedi: {e})")

    print("\n---------------------------------------------------")
    print("Beklenen Çıktı:")
    print("1. Asetik Asit (CC(=O)O)")
    print("2. Salisilik Asit (Oc1ccccc1C(=O)O)")


if __name__ == "__main__":
    main()