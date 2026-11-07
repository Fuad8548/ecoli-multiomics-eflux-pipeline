import os
import cobra

# 1. DEFINE PATHS
CONSTRAINED_MODEL_PATH = "results/models/iJO1366_constrained.xml"

if not os.path.exists(CONSTRAINED_MODEL_PATH):
    raise FileNotFoundError(
        f"Could not find the model file at '{CONSTRAINED_MODEL_PATH}'. "
        "Please ensure Milestone 2 code has executed successfully first."
    )

# 2. LOAD THE PERMANENTLY BURNED E-FLUX MODEL
print(f"Loading constrained SBML model from: {CONSTRAINED_MODEL_PATH}...")
model = cobra.io.read_sbml_model(CONSTRAINED_MODEL_PATH)
print(f"Model successfully charged. Reactions: {len(model.reactions)}, Genes: {len(model.genes)}")

# 3. CONFIGURE ENVIRONMENTAL SUBSTRATE TRAFFIC
# Set up standard physiological glucose uptake rate (-10 mmol/gDCW/h)
# In COBRA, lower_bound < 0 means uptake into the cellular network.
if "EX_glc__D_e" in model.reactions:
    model.reactions.EX_glc__D_e.lower_bound = -10.0
    model.reactions.EX_glc__D_e.upper_bound = 0.0

    if "EX_o2_e" in model.reactions:
        model.reactions.EX_o2_e.lower_bound = -5.0  # Restricts max oxygen intake to 5 mmol/gDCW/h
        model.reactions.EX_o2_e.upper_bound = 0.0
        print(" └── [FIX] Oxygen availability constrained to MICROAEROBIC levels (-5.0).")


# Ensure alternative fermentation secretion sinks are wide open (solver defaults)
# These reactions allow the model to expel accumulated fermentative waste.
fermentation_exchanges = ["EX_ac_e", "EX_lac__D_e", "EX_etoh_e", "EX_for_e", "EX_succ_e"]
for ex_id in fermentation_exchanges:
    if ex_id in model.reactions:
        model.reactions.get_by_id(ex_id).upper_bound = 1000.0  # Open excretion ceiling

# 4. EXECUTE FLUX BALANCE ANALYSIS (FBA)
# iJO1366 uses the Core Biomass Objective Function 'BIOMASS_Ec_iJO1366_core_53_95M'
print(f"\nActive Optimization Target: {model.objective.expression}")
print("Running linear programming optimization...")
solution = model.optimize()

# 5. RENDER PREDICTIVE METABOLIC FLOW SUMMARY
print("\n========================================================")
print("             METABOLIC FLOW PREDICTIONS (FBA)           ")
print("========================================================")
print(f"Optimized Biomass Growth Rate: {solution.objective_value:.4f} hr^-1")

# Extract and map exact exchange fluxes
print("\nPredicted Secretion Profile to Extracellular Matrix:")
fermentation_profile = {
    "EX_ac_e":    "Acetate Secretion     ",
    "EX_lac__D_e": "D-Lactate Secretion   ",
    "EX_etoh_e":   "Ethanol Secretion     ",
    "EX_for_e":    "Formate Secretion     ",
    "EX_succ_e":   "Succinate Secretion   ",
    "EX_co2_e":    "Carbon Dioxide Egress "
}

for rxn_id, descriptive_name in fermentation_profile.items():
    if rxn_id in model.reactions:
        flux_value = solution.fluxes[rxn_id]
        # Only print active fluxes or zero trends clearly
        print(f" └── {descriptive_name} ({rxn_id}): {flux_value:8.4f} mmol/gDCW/h")
        
print("========================================================")
