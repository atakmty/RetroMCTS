#RetroMCTS: AI-Driven Retrosynthesis Planner

RetroMCTS is a Python-based Artificial Intelligence agent designed to solve the retrosynthetic planning problem. Using Monte Carlo Tree Search (MCTS) and chemical rules encoded in SMARTS, the agent recursively decomposes complex target molecules (like drugs) into commercially available building blocks.

This project serves as a Proof of Concept (PoC) for an Artificial Intelligence course, demonstrating the application of Search Algorithms, Logical Agents, and Heuristics in the domain of Computational Chemistry.

🧪 Project Overview

In drug discovery, designing a molecule is only half the battle; knowing how to synthesize it in the lab is equally critical. Retrosynthesis involves working backward from a target molecule to simple starting materials.

RetroMCTS treats this as a search problem:

State Space: Chemical molecules.

Actions: Chemical reactions (e.g., Ester Hydrolysis, Amide Cleavage).

Goal: Reach a state where all fragments exist in a "Buyable Chemicals" database.

Algorithm: MCTS balances exploration (trying new reactions) and exploitation (following promising routes) to find a valid synthesis path efficiently.

🚀 Key Features

Monte Carlo Tree Search (MCTS): Implements Selection, Expansion, Simulation, and Backpropagation phases.

UCB1 Optimization: Uses the Upper Confidence Bound formula to manage the exploration-exploitation trade-off.

Chemical Intelligence: Powered by RDKit for molecular manipulation and reaction handling.

Knowledge Base: Uses SMARTS templates to define chemical rules (e.g., breaking Ester, Amide, and Ether bonds).

Strict Fragment Validation: Implements a logic where a reaction is only valid if all side products are commercially available (AND-logic constraint).

🛠️ Installation & Usage

Prerequisites

Python 3.8+

RDKit

Installation

Install RDKit via pip
pip install rdkit


Running the Agent

Clone the repository.

Open mcts_agent.py.

Set your target molecule (SMILES format) in the __main__ block.

Run the script:

python mcts_agent.py


Example Output

=== SUGGESTED SYNTHESIS ROUTE ===
Target: CC(=O)Nc1ccc(O)cc1 (Paracetamol)
  |
  v (Amide Hydrolysis)
Intermediate: Nc1ccc(O)cc1
Status: SOLVED (Visits: 49)

RESULT: SUCCESS! All fragments (Main + Side Products) are commercially available.


🧠 Algorithmic Logic

The agent operates on a simplified AND-OR Graph Search principle:

Selection: Traverses the tree using UCB1 to find the most promising leaf node.

Expansion: Applies known reaction templates to the molecule.

Constraint: If a reaction produces a side product that is NOT in the buyable database, that branch is pruned immediately (Strict Enforcement).

Simulation: Estimates the "synthesizability" of the new fragments using a heuristic (molecular complexity/atom count).

Backpropagation: Updates the visit count and value estimates up the tree.

⚠️ Limitations & Future Work

This project is an academic prototype and has known limitations:

Knowledge Base: Limited to a small set of reaction rules (Ester, Amide, Ether) to keep the search space manageable without a policy network.

Heuristics: Uses a naive reward function based on atom counts instead of learned Synthetic Accessibility Scores (SAScore).

Reaction Conditions: Ignores physical constraints like temperature, pressure, or stereochemistry.

For a full list of technical constraints and roadmap, please refer to project_limitations.md.

📚 References

Segler, M. H., Preuss, M., & Waller, M. P. (2018). Planning chemical syntheses with deep neural networks and symbolic AI. Nature.

Browne, C. B., et al. (2012). A Survey of Monte Carlo Tree Search Methods. IEEE Transactions on Computational Intelligence and AI in Games.

RDKit: Open-source cheminformatics. https://www.rdkit.org

Author: Ata Kamutay
Course: AIE623 Artificial Intelligence
Date: December 2025
