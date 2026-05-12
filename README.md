# 🧬 Genomic Anomaly Atlas

Welcome to the **Genomic Anomaly Atlas** — a professional analytical pipeline and curated collection of genetic variations. This project bridges the gap between raw DNA sequences and clinical insights using AI-driven automation and visual data analysis.

## 🚀 Key Features: v4.0 "The Visual Architect"
* **Visual Genomic Landscape:** Automated generation of GC-content distribution plots (Radar) for structural analysis.
* **Exon Homology Search:** Comparative engine to identify genetic isoforms (e.g., FGFR2 IIIb vs IIIc).
* **Multi-Level Scanning:** Detection of cancer markers (TERT), neurodegenerative repeats (HTT), and longevity variants (SIRT1).
* **Live API Integration:** Real-time connection to **NCBI ClinVar** for official clinical validation.
* **Professional Reporting:** Automated PDF generation with integrated charts and clinical summaries.

---

## 🏛 Project Structure

### 1. Splicing & Isoforms
* **[FGFR2 Analysis]:** Disruption of epithelial/mesenchymal isoforms.
* **[Exon Similarity Module]:** Automated homology calculation between sample DNA and reference sequences.

### 2. Genetic Anomalies & "Superpowers"
* **[Huntington’s Disease (HTT)]:** CAG trinucleotide repeat expansion analysis.
* **[Sirtuin 1 (SIRT1)]:** Longevity variant (rs7069102) identification.
* **[LRP5 (D171V)]:** The "Titan Bone" mutation (rs121908675) analysis.

---

## 🛠 Installation & Usage

### 1. Requirements
* **Step 1:** Install **Python 3.8+**.
* **Step 2:** Install dependencies:
    ```bash
    pip install biopython reportlab tqdm matplotlib
    ```
* **Step 3:** For correct PDF generation (Cyrillic support), install **DejaVu fonts**:
    ```bash
    sudo apt-get install -y fonts-dejavu-core
    ```

### 2. Quick Start (CLI)
You can run the pipeline with default settings or specify your own files via command line:
* **Run with defaults:**
    ```bash
    python Genomic_Pipeline.py
    ```
* **Custom Analysis:**
    ```bash
    python Genomic_Pipeline.py --input rs_numbers.txt --dna sample_dna.txt --output MyReport.pdf
    ```

### 🐳 Run with Docker (Recommended)
If you have Docker installed, you don't need to install Python or fonts manually:
1. **Build the image:**
    ```bash
    docker build -t genomic-atlas .
    ```
2. **Run the analysis:**
    ```bash
    docker run -v $(pwd):/app genomic-atlas
    ```

### 3. Interactive Demo (Lab Environment)
* **Current Version:** Click to launch the latest analytical environment: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DrabExistence/genomic-anomaly-atlas/blob/main/Genomic_Anomaly_Atlas.ipynb)
* **Legacy Versions:** For project evolution history, see `Demo_Colab(old).ipynb`.

---

## 💻 Tech Stack
* **Language:** Python
* **Bioinformatics:** Biopython (Entrez API), Regular Expressions for motif seeking.
* **Reporting:** ReportLab (PDF), Matplotlib (Data Viz).
* **Environment:** Docker, Google Colab.

---

## 👨‍🔬 About the Author
I am an **AI-Bioinformatics Architect**. I build tools that transform complex genomic data into actionable biological insights.

*Currently exploring the boundaries between human genomics and AI-driven diagnostics.*
