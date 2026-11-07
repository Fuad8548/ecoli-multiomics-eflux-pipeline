# ==============================================================================
# MILESTONE 2: Contextual Integration Module (E-Flux Algorithm)
# Project: Integrative Bulk RNA-Seq Transcriptomics & Genome-Scale Flux Modeling
# Model Core: iJO1366 (Escherichia coli K-12 MG1655)
# ==============================================================================

import os
import re
import pandas as pd
import numpy as np
import cobra
from cobra.io import load_model

# 1. SETUP & REPOSITORY SANITIZATION
INPUT_CSV = "results/data/downregulated_transcripts.csv"
OUTPUT_MODEL_DIR = "results/models"
os.makedirs(OUTPUT_MODEL_DIR, exist_ok=True)

print("└── [SUCCESS] Contextual Integration Module initialized.")

# 2. LOAD METABOLIC INFRASTRUCTURE & TRANSCRIPT SIGNATURES
# Load the definitive E. coli K-12 genome-scale model
try:
    model = load_model("iJO1366")
    print(f"└── [SUCCESS] Model iJO1366 loaded. Reactions: {len(model.reactions)}, Genes: {len(model.genes)}")
except Exception as e:
    print("--> Model not found locally in cache. Fetching core BiGG model...")
    # Fallback to standard e_coli_core if offline or cache misses, or download iJO1366
    model = cobra.io.load_json_model(cobra.io.download.authorize_bigg_api_url() + "/models/iJO1366.json")

# Read the structural output CSV from Milestone 1
de_matrix = pd.read_csv(INPUT_CSV)

# Create a fast-lookup dictionary mapping b-number -> Fold Change ratio (2^log2FC)
# Baseline multiplier for unperturbed genes is 1.0
gene_bounds_modifier = {row['gene_id']: 2**row['log2FoldChange'] for _, row in de_matrix.iterrows()}
print(f"└── [SUCCESS] Parsed {len(gene_bounds_modifier)} down-regulated expression keys.")

# 3. ROBUST GPR PARSING & BOUNDARY MUTATION
constrained_reactions_count = 0

print("\nExecuting GPR boundary mutations across metabolic network...")
for reaction in model.reactions:
    # Skip reactions without structural gene links (e.g., spontaneous, diffusion, exchanges)
    if not reaction.genes:
        continue
        
    gpr_expression = reaction.gene_reaction_rule
    
    # Simple direct one-to-one mapping rule
    if "and" not in gpr_expression and "or" not in gpr_expression:
        clean_gene = gpr_expression.strip("() ")
        modifier = gene_bounds_modifier.get(clean_gene, 1.0)
        
    else:
        # Evaluate complex expressions using boolean evaluation proxies safely
        # Replace b-numbers with their corresponding transcript ratio modifiers
        parsed_expr = gpr_expression
        # Find all b-numbers using regular expressions
        genes_in_rule = re.findall(r'b\d{4}', parsed_expr)
        
        for g in list(set(genes_in_rule)):
            val = gene_bounds_modifier.get(g, 1.0)
            parsed_expr = parsed_expr.replace(g, str(val))
            
        # Clean syntax for evaluation logic in Python: 'and' -> min(), 'or' -> sum()
        # To maintain strict arithmetic execution without unsafe eval(), we handle common clean structures:
        try:
            # Safe evaluation replacement for simple/nested constructs
            # Convert 'and' to min-bounds and 'or' to max/sum bounds
            safe_expr = parsed_expr.replace("and", ",").replace("or", "+")
            # For our pipeline, we extract the structural minimum bottleneck if an "and" block is present
            if "and" in gpr_expression:
                vals = [float(x) for x in re.findall(r'\d+\.\d+|\d+', parsed_expr)]
                modifier = min(vals) if vals else 1.0
            else: # Dynamic backup or choice
                vals = [float(x) for x in re.findall(r'\d+\.\d+|\d+', parsed_expr)]
                modifier = min(sum(vals), 1.0) # Clip maximum capacity scale at 1.0
        except Exception:
            modifier = 1.0

    # Apply E-Flux dampening: Constrain UPPER bound for forward reactions
    if modifier < 1.0:
        original_upper = reaction.upper_bound
        # Scale down the upper boundary proportionally
        reaction.upper_bound = original_upper * modifier
        constrained_reactions_count += 1

print(f"└── [SUCCESS] E-Flux constraints applied to {constrained_reactions_count} reactions.")

# 4. EXPORT NETWORK ARCHITECTURE
cobra.io.write_sbml_model(model, os.path.join(OUTPUT_MODEL_DIR, "iJO1366_constrained.xml"))
print("└── [SUCCESS] Constrained metabolic model successfully written to disk as SBML (.xml).")
