import xml.etree.ElementTree as ET
from urllib.error import URLError, HTTPError
from tqdm.auto import tqdm
from Bio import Entrez
import re
import os
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# --- НАСТРОЙКИ NCBI ---
Entrez.email = "genomic.atlas.project@example.com"

# --- ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ---

def calculate_similarity(seq1, seq2):
    """Рассчитывает процент сходства между двумя последовательностями."""
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    similarity = (matches / max(len(seq1), len(seq2))) * 100
    return round(similarity, 2)

def generate_gc_plot(sequence, window_size=20, output_img="gc_plot.png"):
    """Создает график распределения GC-состава."""
    values = []
    for i in range(len(sequence) - window_size + 1):
        subseq = sequence[i:i + window_size]
        gc = (subseq.count('G') + subseq.count('C')) / len(subseq) * 100
        values.append(gc)

    plt.figure(figsize=(8, 3))
    plt.plot(values, color='#2c3e50', linewidth=1.5, label='GC %')
    plt.fill_between(range(len(values)), values, color='#3498db', alpha=0.3)
    plt.axhline(y=50, color='r', linestyle='--', alpha=0.3)
    plt.title("Genomic Landscape: GC-Content Radar", fontsize=12)
    plt.xlabel("Position (bp)")
    plt.ylabel("GC %")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_img, dpi=100)
    plt.close()
    return output_img

# --- ОСНОВНЫЕ БЛОКИ АНАЛИЗА ---

def get_clinvar_significance_for_list(rs_ids_list):
    """Пакетная проверка rs-номеров через NCBI ClinVar."""
    results = {}
    clinvar_ids = []
    id_map = {}
    if not rs_ids_list: return results

    print(f"🔍 Запрос ClinVar для {len(rs_ids_list)} вариантов...")
    for rs_id in tqdm(rs_ids_list, desc="Поиск в базе"):
        try:
            handle = Entrez.esearch(db="clinvar", term=f"{rs_id}[Variant ID]")
            search_res = Entrez.read(handle)
            handle.close()
            if search_res['IdList']:
                cid = search_res['IdList'][0]
                clinvar_ids.append(cid)
                id_map[cid] = rs_id
            else:
                results[rs_id] = "Not Found"
        except Exception as e:
            results[rs_id] = f"Error"

    if clinvar_ids:
        try:
            handle = Entrez.esummary(db="clinvar", id=",".join(clinvar_ids), validate=False)
            summaries = Entrez.read(handle, validate=False)
            handle.close()
            for doc in summaries['DocumentSummarySet']['DocumentSummary']:
                uid = str(doc.get('uid') or doc.attributes.get('uid'))
                orig_rs = id_map.get(uid)
                if orig_rs:
                    results[orig_rs] = doc['germline_classification']['description']
        except:
            pass
    return results

def universal_genomic_scanner(sequence):
    """Анализ сырой ДНК: GC-состав, повторы и SIRT1."""
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    htt_repeats_count = 0
    htt_repeats_status = "N/A"
    
    repeats = re.findall(r'(?:CAG){3,}', sequence)
    if repeats:
        htt_repeats_count = len(max(repeats, key=len)) // 3
        htt_repeats_status = "PATHOGENIC" if htt_repeats_count >= 36 else "NORMAL"

    sirt1_detected = sequence.startswith('G')
    
    return {
        "gc_content": gc,
        "htt_repeats_count": htt_repeats_count,
        "htt_repeats_status": htt_repeats_status,
        "sirt1_detected": sirt1_detected
    }

def save_report_to_pdf(results, dna_sequence, scanner_data, exon_results, filename="Genomic_Report.pdf"):
    """Генерация финального PDF отчета."""
    try:
        pdfmetrics.registerFont(TTFont('DejaVuSans', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
        font_name = 'DejaVuSans'
    except:
        font_name = 'Helvetica'

    c = canvas.Canvas(filename, pagesize=letter)
    width, height = letter

    # Заголовок
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 18)
    c.drawCentredString(width / 2.0, height - 50, "Genomic Anomaly Atlas: Diagnostic Report")

    # График
    img_path = generate_gc_plot(dna_sequence)
    c.drawImage(img_path, 50, height - 260, width=500, height=160)

    # Результаты сканирования
    c.setFont(font_name, 10)
    curr_y = height - 280
    c.drawString(50, curr_y, f"Global GC-Content: {scanner_data['gc_content']:.2f}%")
    curr_y -= 15
    c.drawString(50, curr_y, f"HTT Gene: {scanner_data['htt_repeats_count']} CAG repeats ({scanner_data['htt_repeats_status']})")
    curr_y -= 15
    c.drawString(50, curr_y, f"SIRT1 Longevity: {'Detected' if scanner_data['sirt1_detected'] else 'Not Detected'}")
    
    # Сравнение экзонов
    curr_y -= 30
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 12)
    c.drawString(50, curr_y, "Exon Comparison Analysis:")
    c.setFont(font_name, 10)
    for name, data in exon_results.items():
        curr_y -= 15
        c.drawString(50, curr_y, f"{name}: {data['similarity']}% Similarity")

    # Таблица ClinVar
    curr_y -= 30
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 12)
    c.drawString(50, curr_y, "ClinVar Clinical Significance:")
    curr_y -= 20
    for rs, status in results.items():
        if rs != "INFO_MESSAGE":
            c.drawString(50, curr_y, f"{rs.upper()}: {status}")
            curr_y -= 15

    c.save()
    print(f"📄 Report saved: {filename}")

# --- ТОЧКА ВХОДА ---

if __name__ == "__main__":
    # Загрузка ДНК
    try:
        with open('sample_dna.txt', 'r') as f: dna = f.read().strip()
    except: dna = "ATGC" * 100
    
    # Сравнение экзонов (пример)
    exon_8 = "ATTATAGTAGAGAGA"
    sim_8 = calculate_similarity(dna[:len(exon_8)], exon_8)
    exon_results = {"Exon_8_IIIb": {"similarity": sim_8}}
    
    # Запуск
    scan_data = universal_genomic_scanner(dna)
    cv_results = get_clinvar_significance_for_list(["rs1137282", "rs121908675"])
    save_report_to_pdf(cv_results, dna, scan_data, exon_results)
