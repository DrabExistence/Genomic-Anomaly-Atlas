# 🧬 Genomic Anomaly Atlas

Welcome to the **Genomic Anomaly Atlas** — a professional analytical pipeline and curated collection of genetic variations. This project bridges the gap between raw DNA sequences and clinical insights using AI-driven automation.

## 🚀 Key Features: v3.0 "The Architect Update"
* **Multi-Level Scanning:** Automated detection of cancer markers (TERT), neurodegenerative repeats (HTT), and "Superhuman" variants (SIRT1, LRP5).
* **Live API Integration:** Real-time connection to **NCBI ClinVar** to fetch official clinical significance.
* **Batch Processing:** Ability to process hundreds of variants from a simple text file.
* **Clinical Reporting:** Automated generation of professional PDF diagnostic reports.

---

## 🏛 Atlas Structure

### 1. Splicing & Point Mutations
* **[FGFR2 Analysis](FGFR2-Analysis.md):** Disruption of epithelial/mesenchymal isoforms.
* **[MTHFR Polymorphism](MTHFR_Analysis.md):** C677T variant and folate cycle impact.

### 2. Genetic Anomalies & "Superpowers"
* **[Huntington’s Disease (HTT)](Huntington_Disease.md):** CAG trinucleotide repeat analysis.
* **[Sirtuin 1 (SIRT1)](Longevity_sirt1.md):** Longevity variants (rs7069102).
* **[LRP5 (D171V)](Titan_bones.md):** The "Titan Bone" mutation (rs121908675) for bone density.

### 3. Oncology & Regulatory Analysis
* **[TERT Promoter](TERT_Promoter_Cancer.md):** Investigating the "immortality switch" in cancer cells.
* **[Visual Genome Analysis](Visual_Genome_Analysis.md):** GC-content distribution and promoter identification.

---

## 🛠 Installation & Usage

### 1. Requirements
* **Step 1:** Install **Python 3.8+**.
* **Step 2:** Ensure you have **External Libraries:** `biopython`, `reportlab`, `tqdm`.
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
* **Step 3:** Add your target RS-numbers to `rs_numbers.txt`.
* **Step 4:** Run the master pipeline:
    ```bash
    python Genomic_Pipeline.py
    ```

### 3. Interactive Demo
* **Step 1:** Open the **[Demo Notebook](https://colab.research.google.com/github/DrabExistence/genomic-anomaly-atlas/blob/main/Demo_Colab.ipynb)** in Google Colab.
* **Step 2:** Follow the cell-by-cell instructions for visualization.
* **Step 3:** Export your results to PDF directly from the browser.

---

## 💻 Tech Stack
* **Language:** Python
* **Bioinformatics:** Biopython (Entrez API), Regular Expressions for motif seeking.
* **Reporting:** ReportLab (PDF Generation).
* **Automation:** Batch processing & API validation.

---

## 👨‍🔬 About the Author
I am an **AI-Bioinformatics Architect**. I build tools that transform complex genomic data into actionable biological insights.

*Currently exploring the boundaries between human genomics and AI-driven diagnostics.*
