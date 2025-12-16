import math
import random
import sys

# RDKit kütüphanesi kontrolü
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("HATA: RDKit kütüphanesi bulunamadı. Lütfen 'pip install rdkit' ile kurun.")
    sys.exit(1)

# --- 1. BİLGİ BANKASI (KNOWLEDGE BASE) ---

REACTION_TEMPLATES = [
    # Kural 1: Ester Hidrolizi (Ester -> Asit + Alkol/Fenol)
    # Örnek: Aspirin'i parçalar.
    {'name': 'Ester Hydrolysis', 'smarts': '[C:1](=[O:2])-[O:3]-[#6:4]>>[C:1](=[O:2])O.[#6:4][O:3]'},

    # Kural 2: Amide Hidrolizi (Amide -> Asit + Amin)
    # Örnek: Parasetamol'ü parçalar.
    # Anlamı: Karbonil(1)-Azot(3) bağını kopar. Karbona OH ekle (Asit), Azotu serbest bırak (Amin).
    {'name': 'Amide Hydrolysis', 'smarts': '[C:1](=[O:2])-[N:3]-[#6:4]>>[C:1](=[O:2])O.[#6:4][N:3]'},

    # Kural 3: Eter Parçalama (Eter -> Alkol + Alkol/Halojenür) - DÜZELTİLDİ
    # Hata Düzeltmesi: Eski kuralda [O:2] sağ tarafta yoktu (unmapped atom).
    # Yeni kural: [O:2] sağ taraftaki ilk parçanın oksijeni olarak atandı ([OH:2]).
    {'name': 'Ether Cleavage', 'smarts': '[#6:1]-[O:2]-[#6:3]>>[#6:1][OH:2].[#6:3]O'}
]

# Hedef (Satın Alınabilir) Moleküller
# Bu listeyi program başladığında normalize edeceğiz.
RAW_BUYABLES = {
    "CC(=O)O",  # Asetik Asit (Aspirin ve Parasetamol parçası)
    "O",  # Su
    "N",  # Amonyak
    "c1ccccc1",  # Benzen
    "CO",  # Metanol
    "CCO",  # Etanol
    "Oc1ccccc1C(=O)O",  # Salisilik Asit (Aspirin parçası)
    "Oc1ccc(N)cc1"  # 4-Aminofenol (Parasetamol parçası - YENİ)
}


# --- 2. YARDIMCI FONKSİYONLAR ---

def canonicalize(smiles):
    """Bir SMILES kodunu RDKit standart formatına çevirir."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    return None


# Satın alınabilir listesini normalize et (Hataları önlemek için)
BUYABLE_CHEMICALS = set()
for s in RAW_BUYABLES:
    norm = canonicalize(s)
    if norm: BUYABLE_CHEMICALS.add(norm)


def is_buyable(smiles):
    """Molekül hedef listede var mı?"""
    norm = canonicalize(smiles)
    return norm in BUYABLE_CHEMICALS


# --- 3. MCTS SINIFLARI ---

class MCTSNode:
    def __init__(self, smiles, parent=None, rule_used=None):
        self.smiles = smiles
        self.parent = parent
        self.rule_used = rule_used

        self.children = []
        self.visits = 0
        self.value = 0.0  # Toplam ödül

        # Durum kontrolü
        self.is_solved = is_buyable(smiles)
        self.is_terminal = self.is_solved  # Eğer çözüldüyse orası bir uçtur

    def ucb_score(self, exploration_weight=1.41):
        """
        Upper Confidence Bound (UCB1) formülü.
        Hem sömürü (exploitation) hem keşif (exploration) dengesini kurar.
        """
        if self.visits == 0:
            return float('inf')  # Hiç ziyaret edilmediyse hemen seç!

        # Ebeveynin ziyaret sayısı (logaritma içinde kullanılır)
        parent_visits = self.parent.visits if self.parent else 1

        exploitation = self.value / self.visits
        exploration = exploration_weight * math.sqrt(math.log(parent_visits) / self.visits)

        return exploitation + exploration

    def expand(self):
        """Olası tüm reaksiyonları uygular ve çocuk düğümleri yaratır."""
        mol = Chem.MolFromSmiles(self.smiles)
        if not mol: return

        # Zaten genişletilmişse tekrar yapma
        if self.children: return

        for action in REACTION_TEMPLATES:
            try:
                rxn = AllChem.ReactionFromSmarts(action['smarts'])
                products_list = rxn.RunReactants((mol,))

                for products in products_list:
                    # Basitleştirme: En büyük parçayı alıyoruz (Ana iskelet)
                    largest_frag = max(products, key=lambda m: m.GetNumAtoms())

                    try:
                        Chem.SanitizeMol(largest_frag)
                        child_smiles = Chem.MolToSmiles(largest_frag, isomericSmiles=True)

                        # Döngüye girmeyi engelle (Anne = Çocuk durumu)
                        if child_smiles == self.smiles: continue

                        child = MCTSNode(smiles=child_smiles, parent=self, rule_used=action['name'])
                        self.children.append(child)
                    except:
                        continue
            except:
                continue


class MCTSEngine:
    def __init__(self, root_smiles):
        self.root = MCTSNode(root_smiles)

    def select(self, node):
        """
        Ağaçta aşağı inerek en iyi (UCB) düğümü seçer.
        Strateji:
        1. Eğer düğümün çocukları yoksa ve genişletilebilir ise -> O düğümü döndür.
        2. Eğer çocukları varsa -> En yüksek UCB puanlı çocuğa git ve tekrar et.
        """
        current = node
        while True:
            # Eğer bu düğüm daha önce hiç genişletilmediyse (ve çözülmediyse), onu seç
            if not current.children and not current.is_terminal:
                # Ancak belki henüz expand edilmemiştir, expand edilebilir mi diye bakmıyoruz
                # MCTS'de genelde: Leaf node'a gel -> Expand et -> Bir çocuğa git.
                return current

            # Eğer terminal ise (çözüldü veya çıkmaz sokak)
            if current.is_terminal or not current.children:
                return current

            # Çocuklardan en yüksek UCB puanına sahip olanı seç
            current = max(current.children, key=lambda c: c.ucb_score())

    def simulate(self, node):
        """
        Rollout: Düğümün ne kadar 'iyi' olduğunu tahmin et.
        Basit bir ödül fonksiyonu:
        - Eğer satın alınabilir ise: +100 puan
        - Değilse: -1 puan (Maliyet)
        """
        if node.is_solved:
            return 100.0
        else:
            # Basit sezgisel: Molekül ne kadar küçükse o kadar iyidir (Basitliğe yaklaşmışızdır)
            # Bu, ajanı molekülü küçültmeye teşvik eder.
            mol = Chem.MolFromSmiles(node.smiles)
            atom_count = mol.GetNumAtoms() if mol else 100
            return 10.0 / atom_count  # Deneysel bir ödül

    def backpropagate(self, node, reward):
        """Sonucu yukarı taşı."""
        current = node
        while current:
            current.visits += 1
            current.value += reward
            current = current.parent

    def run(self, iterations=10):
        print(f"--- MCTS Başlıyor ({iterations} iterasyon) ---")

        for i in range(iterations):
            # 1. Selection
            selected_node = self.select(self.root)

            # 2. Expansion (Eğer seçilen düğüm daha önce ziyaret edildiyse ve terminal değilse genişlet)
            # Not: Standart MCTS'de genelde ilk ziyarette expand edilir.
            if selected_node.visits > 0 and not selected_node.is_terminal:
                selected_node.expand()
                # Genişledikten sonra çocuklardan birini seçmemiz lazım (veya ilkini)
                if selected_node.children:
                    selected_node = selected_node.children[0]  # Basitlik için ilk çocuğa geç

            # 3. Simulation
            reward = self.simulate(selected_node)

            # 4. Backpropagation
            self.backpropagate(selected_node, reward)

            # Log (İsteğe bağlı)
            # print(f"Iter {i+1}: Node={selected_node.smiles[:10]}... Reward={reward:.2f}")

        print("--- Arama Tamamlandı ---")

    def print_best_path(self):
        """En çok ziyaret edilen yolu yazdırır."""
        current = self.root
        print("\n=== ÖNERİLEN SENTEZ YOLU ===")
        print(f"Hedef: {current.smiles}")

        while current.children:
            # En iyi hamle genelde 'en çok ziyaret edilen' (robust) olandır, en yüksek puanlı değil.
            best_child = max(current.children, key=lambda c: c.visits)
            print(f"  |\n  v ({best_child.rule_used})")
            print(f"Ara Ürün: {best_child.smiles} (Ziyaret: {best_child.visits}, Çözüldü: {best_child.is_solved})")

            if best_child.is_solved:
                print("\nSONUÇ: Başarılı! Bu molekül satın alınabilir listesinde.")
                return
            current = best_child

        print("\nSONUÇ: Tam yol bulunamadı veya daha fazla iterasyon gerekli.")


# --- 4. TEST ---
if __name__ == "__main__":
    print("--- MCTS Çoklu Kural Testi ---")

    # Test 1: Aspirin (Ester içerir)
    target_aspirin = "CC(=O)Oc1ccccc1C(=O)O"

    # Test 2: Parasetamol (Amid içerir)
    # Yapı: Asetik asit + 4-Aminofenol birleşimi
    target_paracetamol = "CC(=O)Nc1ccc(O)cc1"

    print(f"\n1. HEDEF: Parasetamol ({target_paracetamol})")
    print("Beklenti: Amid bağının kopması ve 'Asetik Asit' + '4-Aminofenol' çıkması.")

    # Motoru Parasetamol için Başlat
    engine = MCTSEngine(target_paracetamol)

    # Aramayı Başlat
    engine.run(iterations=50)

    # Sonucu Göster
    engine.print_best_path()