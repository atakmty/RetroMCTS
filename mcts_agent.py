import math
import random
import sys
import time
from functools import lru_cache
from collections import Counter

# RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
except ImportError:
    print("ERROR: RDKit library not found.")
    sys.exit(1)

# PubChemPy
try:
    import pubchempy as pcp

    PUBCHEM_AVAILABLE = True
except ImportError:
    PUBCHEM_AVAILABLE = False

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')


# --- COLOR CLASS (For Terminal) ---
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# ==============================================================================
# SECTION 1: CHEMISTRY ENGINE
# ==============================================================================

class ChemistryEngine:
    def __init__(self):
        # Retrosynthetic Rules
        self.rules = [
            {'name': 'Ester Cleavage', 'smarts': '[C:1](=[O:2])-[O:3]-[#6:4]>>[C:1](=[O:2])O.[#6:4][O:3]'},
            {'name': 'Amide Cleavage', 'smarts': '[C:1](=[O:2])-[N:3]-[#6:4]>>[C:1](=[O:2])O.[#6:4][N:3]'},
            {'name': 'Ether Cleavage', 'smarts': '[#6:1]-[O:2]-[#6:3]>>[#6:1][OH:2].[#6:3]O'},
            {'name': 'Imine Cleavage', 'smarts': '[C:1]=[N:2]>>[C:1]=O.[N:2]'},
            {'name': 'Alkene Cleavage', 'smarts': '[C:1]=[C:2]>>[C:1]=O.[C:2]=O'},
            {'name': 'Nitro Reduction (Retro)', 'smarts': '[c:1][NH2]>>[c:1][N+](=O)[O-]'},
            {'name': 'Sulfonamide Cleavage', 'smarts': '[S:1](=[O:2])(=[O:3])-[N:4]>>[S:1](=[O:2])(=[O:3])O.[N:4]'}
        ]

        self.breakable_patterns = [
            Chem.MolFromSmarts('[C](=[O])-[O]-[#6]'),  # Ester
            Chem.MolFromSmarts('[C](=[O])-[N]-[#6]'),  # Amide
            Chem.MolFromSmarts('[#6]-[O]-[#6]'),  # Ether
            Chem.MolFromSmarts('[C]=[N]'),  # Imine
            Chem.MolFromSmarts('[C]=[C]'),  # Alkene
            Chem.MolFromSmarts('[S](=[O])(=[O])-[N]'),  # Sulfonamide
            Chem.MolFromSmarts('[c][NH2]')  # Aromatic Amine
        ]

        self.reactions = []
        for r in self.rules:
            try:
                rxn = AllChem.ReactionFromSmarts(r['smarts'])
                self.reactions.append((r['name'], rxn))
            except:
                pass

    def get_atom_count(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return 9999
        return mol.GetNumAtoms()

    def is_simple_building_block(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return False

        atom_count = mol.GetNumAtoms()
        ring_count = mol.GetRingInfo().NumRings()

        # 1. MICRO MOLECULES
        if atom_count <= 6: return True

        # 2. BREAKABLE BOND CHECK
        has_breakable_bond = False
        for pattern in self.breakable_patterns:
            if mol.HasSubstructMatch(pattern):
                has_breakable_bond = True
                break
        if has_breakable_bond: return False

        # 3. GENERAL PHYSICAL LIMITS
        if ring_count == 1 and atom_count <= 15: return True
        if ring_count == 0 and atom_count <= 15: return True
        if ring_count == 2 and atom_count <= 11: return True

        return False

    def get_valid_moves(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return []
        possible_moves = []

        for name, rxn in self.reactions:
            try:
                outcomes = rxn.RunReactants((mol,))
                for outcome in outcomes:
                    fragments = []
                    valid_outcome = True
                    for frag in outcome:
                        try:
                            Chem.SanitizeMol(frag)
                            smi = Chem.MolToSmiles(frag, isomericSmiles=True)
                            fragments.append(smi)
                        except:
                            valid_outcome = False; break

                    if not valid_outcome or not fragments: continue

                    fragments.sort(key=lambda x: self.get_atom_count(x), reverse=True)
                    main_product = fragments[0]
                    side_products = fragments[1:]

                    if all(self.is_simple_building_block(sp) for sp in side_products):
                        if main_product != smiles:
                            possible_moves.append({
                                'rule': name,
                                'main': main_product,
                                'sides': side_products
                            })
            except:
                continue
        return possible_moves


# ==============================================================================
# SECTION 2: MCTS BRAIN
# ==============================================================================

class Node:
    def __init__(self, state, parent=None, action_used=None, side_products=None):
        self.state = state
        self.parent = parent
        self.action_used = action_used
        self.side_products = side_products if side_products else []
        self.children = []
        self.visits = 0
        self.score = 0.0
        self.untried_moves = None

    def ucb1(self, c=1.41):
        if self.visits == 0: return float('inf')
        return (self.score / self.visits) + c * math.sqrt(math.log(self.parent.visits) / self.visits)


class MCTSSolver:
    def __init__(self, root_smiles, chemistry_engine):
        self.root = Node(root_smiles)
        self.chem = chemistry_engine

    def solve(self, iterations=200):
        # Progress Bar
        print(f"   {Colors.CYAN}Analyzing:{Colors.ENDC} ", end="")
        bar_length = 30

        for i in range(iterations):
            if i % (iterations // bar_length) == 0:
                progress = (i + 1) / iterations
                filled = int(bar_length * progress)
                bar = '█' * filled + '░' * (bar_length - filled)
                sys.stdout.write(f"\r   {Colors.CYAN}Progress:{Colors.ENDC} [{bar}] %{int(progress * 100)}")
                sys.stdout.flush()

            leaf = self.select(self.root)
            child = self.expand(leaf)
            simulation_score = self.simulate(child)
            self.backpropagate(child, simulation_score)

        sys.stdout.write(f"\r   {Colors.CYAN}Progress:{Colors.ENDC} [{'█' * bar_length}] %100     \n")
        print(f"   {Colors.GREEN}✓ Computation Complete.{Colors.ENDC}\n")

    def select(self, node):
        while node.children:
            node = max(node.children, key=lambda n: n.ucb1())
        return node

    def expand(self, node):
        if node.untried_moves is None:
            if self.chem.is_simple_building_block(node.state):
                node.untried_moves = []
            else:
                node.untried_moves = self.chem.get_valid_moves(node.state)

        if not node.untried_moves: return node

        move = node.untried_moves.pop()
        child_node = Node(
            state=move['main'],
            parent=node,
            action_used=move['rule'],
            side_products=move['sides']
        )
        node.children.append(child_node)
        return child_node

    def simulate(self, node):
        atom_count = self.chem.get_atom_count(node.state)
        if self.chem.is_simple_building_block(node.state): return 100.0
        return 20.0 / atom_count

    def backpropagate(self, node, score):
        while node:
            node.visits += 1
            node.score += score
            node = node.parent

    def get_best_trace(self):
        trace = []
        current = self.root
        while True:
            step_info = {
                'smiles': current.state,
                'action': current.action_used,
                'sides': current.side_products,
                'is_solved': self.chem.is_simple_building_block(current.state)
            }
            trace.append(step_info)
            if not current.children: break
            current = max(current.children, key=lambda c: c.visits)
        return trace


# ==============================================================================
# SECTION 3: INTERFACE
# ==============================================================================

@lru_cache(maxsize=1024)
def get_pubchem_name(smiles):
    if not PUBCHEM_AVAILABLE: return smiles
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds and compounds[0].iupac_name:
            return f"{compounds[0].iupac_name}"
        if compounds and compounds[0].synonyms:
            return f"{compounds[0].synonyms[0]}"
    except:
        pass
    return smiles


def print_banner():
    banner = f"""
    {Colors.BOLD}{Colors.CYAN}
    ######################################################
    #                                                    #
    #      R E T R O - M C T S   A I   A G E N T         #
    #                                                    #
    #      AI-Powered Retrosynthesis Tool                #
    #      v2.0                                          #
    #                                                    #
    ######################################################
    {Colors.ENDC}
    """
    print(banner)


def print_instructions():
    print(f"""
    {Colors.BOLD}{Colors.UNDERLINE}HOW IT WORKS?{Colors.ENDC}
    This agent uses Monte Carlo Tree Search (MCTS) to break down the 
    target molecule into simple components. The search consists of 4 stages:

    1. {Colors.BOLD}Selection:{Colors.ENDC} The most promising reaction path is chosen via UCB1.
    2. {Colors.BOLD}Expansion:{Colors.ENDC} Molecule is fragmented using 7 known chemical rules.
    3. {Colors.BOLD}Simulation:{Colors.ENDC} Fragment complexity is scored based on atom count.
    4. {Colors.BOLD}Backpropagation:{Colors.ENDC} The score is carried back up the tree.

    {Colors.BOLD}{Colors.UNDERLINE}AGENT REACTION RULES (SCISSORS):{Colors.ENDC}
    ✂️  {Colors.CYAN}Ester & Amide Cleavage:{Colors.ENDC} Breaks basic pharmaceutical bonds.
    ✂️  {Colors.CYAN}Ether & Imine Cleavage:{Colors.ENDC} Fragments Oxygen and Nitrogen bridges.
    ✂️  {Colors.CYAN}Alkene Cleavage:{Colors.ENDC} Oxidizes Carbon-Carbon double bonds.
    ✂️  {Colors.CYAN}Nitro Reduction (Retro):{Colors.ENDC} Converts Amines back to Nitro groups.
    ✂️  {Colors.CYAN}Sulfonamide Cleavage:{Colors.ENDC} Breaks down Sulfate-class drugs.
    """)


def main():
    chem = ChemistryEngine()
    print_banner()
    print_instructions()

    if not PUBCHEM_AVAILABLE:
        print(f"{Colors.WARNING}WARNING: 'pubchempy' not found. Compound naming disabled.{Colors.ENDC}")

    while True:
        print(f"{Colors.BOLD}Ready for New Analysis.{Colors.ENDC} (Exit: 'q')")
        smi = input(f">> Target SMILES Code: ").strip()

        if smi.lower() == 'q':
            print("\nGoodbye! 👋")
            break
        if not smi: continue

        # Validation
        mol_check = Chem.MolFromSmiles(smi)
        if mol_check is None:
            print(f"{Colors.FAIL}❌ ERROR: Invalid SMILES format!{Colors.ENDC}\n")
            continue

        try:
            target_name = get_pubchem_name(smi)
            print(f"\n{Colors.BOLD}Target Molecule:{Colors.ENDC} {target_name}")
            print(f"{Colors.BOLD}Formula:{Colors.ENDC} {smi}")
            print("-" * 50)

            solver = MCTSSolver(smi, chem)
            solver.solve(iterations=200)

            trace = solver.get_best_trace()

            # --- REPORTING SECTION ---
            print(f"\n{Colors.BOLD}{Colors.UNDERLINE}SYNTHETIC ROUTE REPORT{Colors.ENDC}\n")

            total_materials = []

            # Root Node Check
            if len(trace) == 1 and trace[0]['is_solved']:
                print(f"{Colors.GREEN}ℹ️  RESULT: This molecule is already a 'Building Block'.{Colors.ENDC}")
                print(f"   (Available in market or too simple to fragment.)\n")
                continue

            # Step-by-Step Printing
            for i, step in enumerate(trace):
                name = get_pubchem_name(step['smiles'])

                if i == 0:
                    # Starting point
                    print(f"🎯 {Colors.BOLD}TARGET:{Colors.ENDC} {name}")
                else:
                    # Arrow and Rule
                    rule = step['action']
                    print(f"   │")
                    print(f"   ▼  {Colors.WARNING}Reaction: {rule}{Colors.ENDC}")

                    # Products
                    if step['is_solved']:
                        color = Colors.GREEN
                        tag = "✅ [BUILDING BLOCK]"
                    else:
                        color = Colors.CYAN
                        tag = "🔹 [INTERMEDIATE]"

                    print(f"   │")
                    print(f"   └── {color}{tag} {name}{Colors.ENDC}")

                    # Side Products
                    if step['sides']:
                        for side in step['sides']:
                            side_name = get_pubchem_name(side)
                            print(f"      └── ➕ {Colors.GREEN}[SIDE PRODUCT] {side_name}{Colors.ENDC}")
                            total_materials.append(side_name)

                # Add final piece to materials if solved
                if i == len(trace) - 1 and step['is_solved']:
                    total_materials.append(name)

            # Final Status
            final_step = trace[-1]
            print("\n" + "=" * 50)

            if final_step['is_solved']:
                print(f"{Colors.GREEN}{Colors.BOLD}✅ SYNTHESIS SUCCESSFULLY PLANNED!{Colors.ENDC}")
                print("All components reduced to basic building blocks.")

                print(f"\n{Colors.BOLD}📃 REQUIRED MATERIALS LIST (Recipe):{Colors.ENDC}")
                counts = Counter(total_materials)
                for mat, count in counts.items():
                    print(f"   • {count}x {mat}")
            elif len(trace) > 1:
                print(f"{Colors.WARNING}{Colors.BOLD}⚠️  PARTIALLY SUCCESSFUL{Colors.ENDC}")
                print("Molecule simplified but could not reach full building blocks.")
                print("Current rules cannot fragment the remaining part.")
            else:
                print(f"{Colors.FAIL}{Colors.BOLD}❌ NO SOLUTION FOUND{Colors.ENDC}")
                print("No suitable starting move found for this molecule.")

            print("=" * 50 + "\n")

        except Exception as e:
            print(f"{Colors.FAIL}An unexpected error occurred: {e}{Colors.ENDC}\n")


if __name__ == "__main__":
    main()