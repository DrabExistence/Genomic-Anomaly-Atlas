import re
import os
import argparse
import matplotlib.pyplot as plt
from Bio import Entrez
from tqdm.auto import tqdm
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# --- НАСТРОЙКИ ---
Entrez.email = "genomic.atlas.project@example.com"

# Константные эталоны из твоего исследования
EXON_8_IIIb = "ATTATAGTAGAGAGA"
EXON_9_IIIc = "GAACTGTGGTATAGA"

# --- МАТЕМАТИКА ---

def calculate_similarity(seq1, seq2):
    """Сравнение последовательностей в процентах."""
    min_len = min(len(seq1), len(seq2))
    if max(len(seq1), len(seq2)) == 0: return 0
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    return round((matches / max(len(seq1), len(seq2))) * 100, 2)

def generate_gc_plot(sequence, window_size=20, output_img="gc_plot.png"):
    """Визуализация ландшафта ДНК."""
    values = [((s.count('G') + s.count('C')) / len(s) * 100) 
              for i in range(len(sequence) - window_size + 1) 
              if (s := sequence[i:i + window_size])]
    
    plt.figure(figsize=(8, 3))
    plt.plot(values, color='#2c3e50', linewidth=1.5)
    plt.fill_between(range(len(values)), values, color='#3498db', alpha=0.3)
    plt.title("Genomic Landscape: GC-Content Radar", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_img, dpi=100)
    plt.close()
    return output_img

# --- АНАЛИЗ ---

def get_clinvar_data(rs_list):
    """Сбор данных из базы ClinVar."""
    results = {}
    if not rs_list: return results
    print(f"🔍 Запрос ClinVar для {len(rs_list)} вариантов...")
    for rs in tqdm(rs_list, desc="NCBI Search"):
        try:
            h = Entrez.esearch(db="clinvar", term=f"{rs}[Variant ID]")
            r = Entrez.read(h); h.close()
            if r['IdList']:
                sum_h = Entrez.esummary(db="clinvar", id=r['IdList'][0])
                sum_r = Entrez.read(sum_h); sum_h.close()
                results[rs] = sum_r['DocumentSummarySet']['DocumentSummary'][0]['germline_classification']['description']
            else: results[rs] = "Not Found"
        except: results[rs] = "API Error"
    return results

# --- ОТЧЕТ ---

def create_pdf(filename, dna, cv_results, scan_data, sim_results):
    try:
        pdfmetrics.registerFont(TTFont('DejaVuSans', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
        f_name = 'DejaVuSans'
    except: f_name = 'Helvetica'

    c = canvas.Canvas(filename, pagesize=letter)
    c.setFont(f"{f_name}-Bold", 16)
    c.drawCentredString(300, 750, "Genomic Anomaly Atlas Report")
    
    # График
    img = generate_gc_plot(dna)
    c.drawImage(img, 50, 530, width=500, height=150)
    
    # Данные
    curr_y = 510
    c.setFont(f"{f_name}-Bold", 12)
    c.drawString(50, curr_y, "I. Primary Metrics:")
    c.setFont(f_name, 10)
    curr_y -= 20
    c.drawString(60, curr_y, f"• GC-Content: {scan_data['gc']:.2f}%")
    c.drawString(250, curr_y, f"• HTT Repeats: {scan_data['htt']} ({scan_data['htt_status']})")
    
    curr_y -= 40
    c.setFont(f"{f_name}-Bold", 12)
    c.drawString(50, curr_y, "II. Exon Similarity (Isoform Analysis):")
    c.setFont(f_name, 10)
    for name, val in sim_results.items():
        curr_y -= 20
        c.drawString(60, curr_y, f"• {name}: {val}% similarity")

    curr_y -= 40
    c.setFont(f"{f_name}-Bold", 12)
    c.drawString(50, curr_y, "III. ClinVar Results:")
    c.setFont(f_name, 9)
    for rs, st in cv_results.items():
        curr_y -= 15
        c.drawString(60, curr_y, f"{rs.upper()}: {st}")
        if curr_y < 50: c.showPage(); curr_y = 750

    c.save()

# --- CLI ---

def main():
    parser = argparse.ArgumentParser(description="Genomic Atlas CLI")
    parser.add_argument("-i", "--input", default="rs_numbers.txt")
    parser.add_argument("-d", "--dna", default="sample_dna.txt")
    parser.add_argument("-o", "--output", default="Genomic_Report.pdf")
    args = parser.parse_args()

    # Загрузка ДНК
    if os.path.exists(args.dna):
        with open(args.dna, 'r') as f: dna = f.read().strip()
    else: dna = "ATGC" * 100 # Default

    # Загрузка RS
    if os.path.exists(args.input):
        with open(args.input, 'r') as f: rs = [l.strip() for l in f if l.strip()]
    else: rs = ["rs1137282", "rs121908675"]

    # Анализ
    scan = {
        "gc": (dna.count('G') + dna.count('C')) / len(dna) * 100,
        "htt": len(max(re.findall(r'(?:CAG){3,}', dna), key=len, default="")) // 3,
        "htt_status": "NORMAL"
    }
    scan["htt_status"] = "PATHOGENIC" if scan["htt"] >= 36 else "NORMAL"
    
    sims = {
        "FGFR2_Exon8_IIIb": calculate_similarity(dna[:len(EXON_8_IIIb)], EXON_8_IIIb),
        "FGFR2_Exon9_IIIc": calculate_similarity(dna[:len(EXON_9_IIIc)], EXON_9_IIIc)
    }
    
    cv = get_clinvar_data(rs)
    create_pdf(args.output, dna, cv, scan, sims)
    print(f"✅ Report saved to {args.output}")

if __name__ == "__main__":
    main()
