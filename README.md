# 🧬 Genomic Anomaly Atlas

Welcome to the **Genomic Anomaly Atlas** — a professional analytical pipeline and curated collection of genetic variations. This project bridges the gap between raw DNA sequences and clinical insights using AI-driven automation.

## 🚀 Key Features: v4.0 "The Visual Architect"
* **Visual Genomic Landscape:** Automated generation of GC-content distribution plots for structural analysis.
* **Exon Similarity Analysis:** Comparative engine to identify genetic isoforms and splicing patterns.
* **Multi-Level Scanning:** Detection of cancer markers (TERT), neurodegenerative repeats (HTT), and longevity variants (SIRT1).
* **Live API Integration:** Real-time connection to **NCBI ClinVar** for clinical validation.
* **Automated Reporting:** Professional PDF generation with integrated data visualization.

---

## 🏛 Atlas Structure

### 1. Splicing & Isoforms
* **[FGFR2 Analysis](FGFR2-Analysis.md):** The switch between epithelial (IIIb) and mesenchymal (IIIc) isoforms.
* **[Similarity Module]:** Automated homology calculation between sample DNA and reference exons.

### 2. Genetic Anomalies & "Superpowers"
* **[Huntington’s Disease (HTT)](Huntington_Disease.md):** CAG trinucleotide repeat expansion analysis.
* **[Longevity & Density]:** Analysis of SIRT1 (rs7069102) and LRP5 (rs121908675) variants.

---

## 🛠 Installation & Usage

### 1. Requirements
* **Step 1:** Install **Python 3.8+**.
* **Step 2:** Ensure you have the necessary **External Libraries**: `biopython`, `reportlab`, `tqdm`, `matplotlib`.
* **Step 3:** For correct PDF generation (Cyrillic support), install **DejaVu fonts**:
    ```bash
    sudo apt-get install -y fonts-dejavu-core
    ```

### 2. Quick Start
* **Step 1:** Clone the repository:
    ```bash
    git clone [https://github.com/DrabExistence/genomic-anomaly-atlas.git](https://github.com/DrabExistence/genomic-anomaly-atlas.git)
    ```
* **Step 2:** Install dependencies using:
    ```bash
    pip install -r requirements.txt
    ```
* **Step 3:** Add your target RS-numbers to `rs_numbers.txt` and raw DNA to `sample_dna.txt`.
* **Step 4:** Run the master pipeline:
    ```bash
    python Genomic_Pipeline.py
    ```
### 🐳 Run with Docker (Recommended)
If you have Docker installed, you don't need to install Python or fonts manually.

* **Step 1: Build the image**
    ```bash
    docker build -t genomic-atlas .
    ```
* **Step 2: Run the analysis**
    ```bash
    docker run -v $(pwd):/app genomic-atlas
    ```
    *Note: The output report (`Genomic_Report.pdf`) will appear in your current folder.*

### 3. Interactive Demo (Lab Environment)
* **Current Version:** Click to launch the latest analytical environment: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DrabExistence/genomic-anomaly-atlas/blob/main/Demo_Colab(new).ipynb)
* **Legacy Versions:** For project evolution history, see `Demo_Colab(old).ipynb`.
* **Instructions:** Follow the cell-by-cell execution in Colab to generate plots and export the final PDF report.

---

## 💻 Tech Stack
* **Language:** Python
* **Bioinformatics:** Biopython (Entrez API), Regular Expressions for motif seeking.
* **Reporting:** ReportLab (PDF Generation), Matplotlib (Data Viz).
* **Automation:** Batch processing & API validation.

---

## 👨‍🔬 About the Author
I am an **AI-Bioinformatics Architect**. I build tools that transform complex genomic data into actionable biological insights.

*Currently exploring the boundaries between human genomics and AI-driven diagnostics.*
