# 1. Code Review & Structural Analysis

## Codebase Walkthrough: Step-by-Step Explanation

The script implements the **E-Flux** algorithm, a well-established method in systems biology used to integrate transcriptomic data directly into a genome-scale metabolic model (GEM). Instead of using hard arbitrary thresholds to classify genes as simply "on" or "off", E-Flux uses the continuous expression levels (fold-changes) as a proxy for maximum enzyme capacity, scaling down the upper flux bounds ($\(V_{max}\)$) of reactions proportionally.
Here is the step-by-step breakdown of how the codebase accomplishes this:

### Step 1: Load Metabolic Infrastructure & Transcript Signatures

```python
try:
    model = load_model("iJO1366")
    ...
except Exception as e:
    model = cobra.io.load_json_model(cobra.io.download.authorize_bigg_api_url() + "/models/iJO1366.json")

de_matrix = pd.read_csv(INPUT_CSV)
gene_bounds_modifier = {row['gene_id']: 2**row['log2FoldChange'] for _, row in de_matrix.iterrows()}
```

- Purpose: Imports the foundational structural network of Escherichia coli metabolism (iJO1366) and map out the experimental data.
- Mechanism:
	1. It attempts to load iJO1366 from a local cache. If it misses, it uses the BiGG API to download the complete model dynamically over JSON.
	2. It reads the differentially expressed (DE) gene dataset using pandas.
	3. It populates a Python dictionary (`gene_bounds_modifier`) by converting the $\(\log_2(\text{Fold Change})\)$ values back into raw linear scale values ($\(2^{\log _{2}FC}\)$). For example, a $\(\log _{2}FC\)$ of -1 converts to a modifier of 0.5 (meaning expression dropped by half).

### Step 2: Robust GPR Parsing & Boundary Mutation (The Core Algorithm)

This block loops over every single reaction in the model to parse its **Gene-Protein-Reaction (GPR)** rules and scale metabolic boundaries.

```python
for reaction in model.reactions:
    if not reaction.genes:
        continue
```
- Reactions mediated by non-protein mechanisms (spontaneous reactions, simple diffusion, or boundary exchange conditions) are ignored since they are independent of transcript levels.

#### Case A: Simple 1-to-1 Mapping

```python
if "and" not in gpr_expression and "or" not in gpr_expression:
    clean_gene = gpr_expression.strip("() ")
    modifier = gene_bounds_modifier.get(clean_gene, 1.0)
```
- If a single reaction is carried out by a single gene product, its scaling factor (modifier) is directly pulled from our expression dictionary. If the gene is not in the down-regulated list, its baseline default remains unperturbed (1.0).

#### Case B: Complex Evaluation (Isoenzymes and Enzyme Complexes)
When multiple genes control a single reaction, the script parses boolean relationships:

```python
genes_in_rule = re.findall(r'b\d{4}', parsed_expr)

for g in list(set(genes_in_rule)):
    val = gene_bounds_modifier.get(g, 1.0)
    parsed_expr = parsed_expr.replace(g, str(val))
```
- Regex Lookup: It looks for systemic E. coli locus tags (e.g., b0001) and swaps out the text string name for its numerical expression ratio modifier calculated in Step 2.

```python
if "and" in gpr_expression:
    vals = [float(x) for x in re.findall(r'\d+\.\d+|\d+', parsed_expr)]
    modifier = min(vals) if vals else 1.0
else: 
    vals = [float(x) for x in re.findall(r'\d+\.\d+|\d+', parsed_expr)]
    modifier = min(sum(vals), 1.0)
```
- Biological Logic Constraints:
	- and (Enzyme Complexes): Subunits working together are limited by a bottleneck principle. If a reaction requires protein A and protein B, the script calculates a `min()` function. The most down-regulated subunit determines the constraint for the whole assembly.
	- or (Isoenzymes): Redundant catalysts working independently pool their capacities together (sum(vals)), reflecting that both enzymes can handle the flux load. The script applies a safety limit (`min(..., 1.0)`) to ensure unperturbed alternative paths don't artifactually balloon the capacity beyond basal limits.

**Constraining the Upper Bounds**
```python
if modifier < 1.0:
    original_upper = reaction.upper_bound
    reaction.upper_bound = original_upper * modifier
    constrained_reactions_count += 1
```
- Wherever a gene is down-regulated (modifier < 1.0), its corresponding reaction maximum speed limit (`upper_bound`) is throttled down proportionally ($\(V_{max, new} = V_{max, old} \times \text{modifier}\)$).

### Step 3: Export Constrained Network Architecture

```python
cobra.io.write_sbml_model(model, os.path.join(OUTPUT_MODEL_DIR, "iJO1366_constrained.xml"))
```
- Purpose: Permanently serializes the newly contextualized, condition-specific metabolic model.
- Mechanism: Saves the updated model object back to an advanced SBML format (`iJO1366_constrained.xml`) with the expression constraints baked into the reaction bounds, making it fully ready for downstream optimization studies like Flux Balance Analysis (FBA).


The script is cleanly written and effectively transitions our pipeline from static data analytics (`DESeq2`) into predictive systems biology (`COBRA`). Here is what it does exceptionally well, along with one crucial mathematical caveat:
- Direct Scaling Layer: It converts the R script's logarithmic data (`log2FoldChange`) back into a natural linear multiplier using 2**row['log2FoldChange']. This changes standard statistical values into a direct biological capacity metric.
- GPR Boolean Logic Integration: It handles complex, multi-subunit enzymes (like `sdhCDAB` or `sucABCD`) using regular expressions (`re.findall`). By applying a `min()` logic pool to "and" associations, it correctly assumes that the most deeply suppressed subunit acts as the structural bottleneck for the entire enzyme complex.
- Proportional Dampening: Rather than brutally knocking out a reaction to zero (which caused the 0.0000 mathematical failure in our previous script), it gracefully multiplies the reaction's maximum theoretical ceiling (`upper_bound`) by the expression ratio (`original_upper * modifier`). This maintains model feasibility.


## 2. The XML output
The script outputs a customized SBML (Systems Biology Markup Language) file (`.xml`). This format is the global standard for computational biology. Its utility in our project includes:
1. Permanence: Instead of re-running the heavy GPR-parsing regex loops every single time we want to test a hypothesis, our dataset constraints are now permanently burned directly into the model's structural file.
2. Universal Compatibility: This specific `.xml` file can now be opened, solved, and visualized across completely different programming environments and software suites. We can pass this file directly into:
	- Python (cobrapy)
	- MATLAB (COBRA Toolbox)
	- Standalone software platforms like Escher (for carbon pathway layout visualization) or OptFlux.

### What does the iJO1366 file structure mean?
Looking inside that XML block explains how the code interacts with the biology:
- The Model Metadata Header:
```xml
<model metaid="meta_iJO1366" id="iJO1366" name="Escherichia coli str. K-12 substr. MG1655" ...>
```
This establishes the structural baseline for E. coli K-12. Our code updates this structural matrix by adding an unlisted constraint layer.

- The Compartment Maps:
```xml
<compartment id="c" name="cytosol" .../>
<compartment id="p" name="periplasm" .../>
```
Unlike simplified core models, iJO1366 models cellular geography explicitly. Rerouting our flux into fermentation will cause protons (H⁺) and products to cross from the cytosol (c), through the periplasm (p), and exit into the extracellular space (e).

- The Metabolite ID Registries (<listOfSpecies>):
```xml
<species metaid="meta_M_akg_c" id="M_akg_c" name="2-Oxoglutarate" ... chemicalFormula="C5H4O5">
```
This is the exact mathematical tracking system. Alpha-Ketoglutarate (akg_c) registered here with its explicit formula. When our code throttles AKGD (Alpha-ketoglutarate dehydrogenase), the model's internal stoichiometric matrix knows exactly how much carbon (C5) is backing up inside the cytoplasm, forcing it into alternative loops.















