# Tetradyme: Tetra-Primer ARMS-PCR Design Suite

**Tetradyme** is a high-performance, web-based bioinformatics tool engineered for the automated thermodynamic design and optimization of inner primers for Tetra-Primer Amplification Refractory Mutation System Polymerase Chain Reaction (T-ARMS-PCR).

---

## Scientific Background

### Principle of Tetra-Primer ARMS-PCR
Tetra-primer ARMS-PCR is a rapid, cost-effective, and highly specific method for genotyping single nucleotide polymorphisms (SNPs). It employs four uniquely designed primers in a single multiplex PCR reaction tube to simultaneously detect both wild-type and mutant alleles.

The system utilizes two primer pairs:
1. **Outer Primers (OF & OR)**: Flank the SNP region, serving as an internal reaction control and generating a large "outer" amplicon.
2. **Inner Primers (IF & IR)**: Designed in opposite orientations. The 3' terminal nucleotide of each inner primer is positioned exactly over the SNP site, making them strictly allele-specific. 

```text
       Outer Forward (OF) ---------------------------> 
5' -----------------------------------------*----------------------------------------- 3'
3' -----------------------------------------*----------------------------------------- 5'
                        <------------------- Inner Reverse (IR) (Allele 1 Specific)
                                            Inner Forward (IF) ----------------------> (Allele 2 Specific)
                                            <------------------------- Outer Reverse (OR)
```
*(Where `*` represents the targeted Single Nucleotide Polymorphism)*

### Mechanism of Action
The technique relies on the principle that Taq DNA polymerase lacks 3' → 5' exonuclease proofreading activity. If the 3' end of an inner primer mismatches the template DNA (i.e., the wrong allele is present), polymerase extension is structurally hindered (refractory). Amplification only proceeds efficiently if the 3' terminal base perfectly matches the target allele.

The outer and inner primers are intentionally engineered asymmetrically around the SNP. This asymmetry guarantees that the two allele-specific amplicons (OF + IR vs. IF + OR) have distinctly different molecular weights, allowing for unambiguous visual differentiation via standard agarose gel electrophoresis.

---

## Application Architecture

Tetradyme features a robust decoupled architecture:
- **Backend (Core Logic)**: Driven by `FastAPI`, `BioPython`, and `primer3-py`. Performs exhaustive thermodynamic simulations, nearest-neighbor melting temperature (Tm) modeling, and secondary structure (hairpin/homodimer) predictions.
- **Frontend (UI)**: A decoupled, responsive Vanilla JS/CSS client that utilizes native `plotly.js` to simulate gel electrophoresis banding patterns based on computed amplicon sizes.

### Algorithmic Enhancements
- **Dynamic Search Space**: The search algorithm evaluates primer sequences ranging from 15 to 35 base pairs, vastly outperforming legacy algorithms restricted to 18-25 bp constraints.
- **Thermodynamic Fallback Logic**: If strict user-defined Tm tolerances cannot be met due to localized AT/GC biases, the algorithm utilizes a greedy heuristic to identify and return the nearest thermodynamically stable primer pair within an expanded ±2.0°C deviation threshold, ensuring a successful design yield for challenging genomic targets.

---

## Deployment & Setup

The repository is fully Dockerized and container-ready for cloud deployments.

### Local Development
To run the server locally, install the core dependencies and launch the backend using Uvicorn:

```bash
# 1. Install pip dependencies
pip install -r requirements.txt

# 2. Start the FastAPI application
python -m uvicorn main:app
```
Once the application startup is complete, the application will be accessible locally via your web browser.

### Cloud Deployment (Render / Hugging Face)
Tetradyme can be deployed as a public web service in minutes:
1. Connect this repository to your preferred hosting provider (e.g., Render, Railway, or Hugging Face Spaces).
2. Select **Docker** as the deployment environment.
3. The provider will automatically parse the included `Dockerfile`, resolve dependencies, and generate your live public URL. 

---

## Data Input Constraints
The processing engine strictly requires two inputs:
1. **Genomic Template**: A standard `.fasta` or `.fa` file containing the target sequences.
2. **SNP Index**: A `.csv` containing exactly five required headers: `ID, pos, ref, alt, outer_tm` where `pos` represents the 1-based index of the polymorphic site.
