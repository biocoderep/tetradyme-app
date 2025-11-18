# Tetradyme - Tetra-Primer ARMS-PCR Designer

Tetradyme is an interactive bioinformatics web application designed to automate and simplify the design of inner primers for Tetra-Primer Amplification Refractory Mutation System PCR (ARMS-PCR). It ensures robust primer design with allele-specific discrimination for SNP genotyping.

## Live Application

The application is deployed and accessible at: https://tinyurl.com/tetradyme

## Overview

Tetra-Primer ARMS-PCR is a widely used molecular technique for SNP genotyping that requires careful primer design to ensure allele-specific amplification. Tetradyme streamlines this process by automating primer design with comprehensive quality control checks, reducing the time and expertise required for manual primer design.

## Understanding Tetra-Primer ARMS-PCR

Tetra-Primer ARMS-PCR is a single-tube PCR method for SNP genotyping that uses four primers in one reaction:

### Primer Types and Functions

**Outer Primers (Forward and Reverse)**
- Flank the SNP region on both sides
- Amplify a control product that confirms successful PCR amplification
- Serve as common primers for both alleles
- Generate the largest amplicon in the reaction

**Inner Primers (Allele-Specific)**
- Two primers designed to discriminate between the two alleles of a SNP
- Each inner primer is specific to one allele (wild-type or mutant)
- Positioned with their 3' ends at or near the SNP site
- The 3' terminal nucleotide matches one of the two SNP alleles
- Work in combination with the opposite outer primer to generate allele-specific products

### PCR Products
The reaction produces three possible outcomes:
- Homozygous wild-type: Control band + Wild-type allele-specific band
- Homozygous mutant: Control band + Mutant allele-specific band  
- Heterozygous: Control band + Both allele-specific bands

## How Tetradyme Designs Inner Primers

Tetradyme focuses specifically on the design of allele-specific inner primers based on SNP location:

### Design Algorithm

1. **SNP Localization**: The application identifies the exact position of the SNP within the uploaded sequence

2. **3' End Positioning**: Inner primers are designed with their 3' terminal nucleotide positioned at the SNP site to ensure maximum allele discrimination

3. **Allele-Specific Design**: 
   - One inner primer is designed with the 3' end matching the reference/wild-type allele
   - The second inner primer is designed with the 3' end matching the alternate/mutant allele
   
4. **Orientation Selection**: 
   - Inner forward primer is paired with the outer reverse primer
   - Inner reverse primer is paired with the outer forward primer

5. **Quality Optimization**: Each primer is evaluated for:
   - Appropriate melting temperature (Tm)
   - GC content within acceptable ranges
   - Presence of GC clamp at the 3' end
   - Absence of secondary structures and primer-dimers
   - Optimal primer length

This design strategy ensures that each inner primer will only efficiently amplify when its 3' end perfectly matches the target allele, enabling reliable SNP discrimination.

## Key Features

- **Streamlined Input Processing**: Accepts FASTA format sequences and CSV-formatted SNP data
- **Automated Primer Design**: Generates allele-specific inner primers optimized for ARMS-PCR
- **Comprehensive Quality Control**:
  - Melting temperature (Tm) calculation using Wallace rule and IDT-compatible settings
  - SNP-specific 3' end positioning for allele discrimination
  - GC clamp validation for primer stability
  - Complete primer property analysis including GC content, molecular weight, and extinction coefficient
- **Interactive Interface**: Built with Streamlit for intuitive user experience
- **Exportable Results**: Download primer designs as CSV files for documentation and further analysis
- **Advanced Validation Options**: Optional secondary structure analysis and primer-dimer detection with NCBI BLAST integration

## Input Requirements

### FASTA File
- Standard FASTA format (.fasta or .fa)
- Contains the target DNA sequence encompassing the SNP region
- Maximum file size: 2000KB

### SNP Data File
- CSV format with required columns (specified in the application interface)
- Contains Single Nucleotide Polymorphism information
- Maximum file size: 2000KB

## Installation and Local Setup

### Prerequisites
- Python 3.7 or higher
- pip package manager

### Installation Steps

1. Clone the repository:
```bash
git clone https://github.com/biocoderep/tetradyme-app.git
cd tetradyme-app
```

2. Create and activate a virtual environment (recommended):
```bash
# Windows
python -m venv venv
venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

3. Install required dependencies:
```bash
pip install -r requirements.txt
```

4. Launch the application:
```bash
streamlit run app.py


## Usage Instructions

1. Access the application through the deployed URL or local server
2. Upload your FASTA sequence file containing the target region
3. Upload the CSV file with SNP information
4. Configure reaction conditions and primer parameters according to your experimental requirements
5. Review the automatically generated primer designs
6. Export results to CSV format for record-keeping and laboratory use

## Authors

**Pooja B.N.** - Primary Developer and Maintainer
GitHub: @biocoderep

**D. Swetha** - Contributing Developer
GitHub: @swetha-pixel


## Citation

If you use Tetradyme in your research or projects, please cite:

```
Pooja B.N. and D.Swetha (2025). Tetradyme: A Web-Based Tool for Tetra-Primer ARMS-PCR Design.
GitHub repository: https://github.com/biocoderep/tetradyme-app
```

## License

This project is licensed under the MIT License. See the LICENSE file in the repository for full details.

## Repository

https://github.com/biocoderep/tetradyme-app

## Contact

For questions, issues, or contributions, please open an issue on the GitHub repository.

---

Tetradyme - Simplifying SNP Genotyping Primer Design
