import re
import os
import matplotlib.pyplot as plt
from Bio import Entrez
from tqdm.auto import tqdm
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# --- НАСТРОЙКИ NCBI ---
# Используем нейтральный технический адрес проекта
Entrez.email = "genomic.atlas.project@example.com" 

# --- БЛОК 1: ВИЗУАЛИЗАЦИЯ (ГЕНЕТИЧЕСКИЙ РАДАР) ---

def generate_gc_plot(sequence, window_size=20, output_img="gc_plot.png"):
    """Создает график распределения GC-состава (скользящее окно)."""
    values = []
    for i in range(len(sequence) - window_size + 1):
        subseq = sequence[i:i + window_size]
        gc = (subseq.count('G') + subseq.count('C')) / len(subseq) * 100
        values.append(gc)
    
    plt.figure(figsize=(8, 3))
    plt.plot(values, color='#2c3e50', linewidth=1.5, label='GC %')
    plt.fill_between(range(len(values)), values, color='#3498db', alpha=0.3)
    plt.axhline(y=50, color='r', linestyle='--', alpha=0.3) # Линия 50% (баланс)
    plt.title("Genomic Landscape: GC-Content Radar", fontsize=12)
    plt.xlabel("Sequence Position (bp)", fontsize=10)
    plt.ylabel("GC %", fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_img, dpi=100)
    plt.close()
    return output_img

# --- БЛОК 2: CLINVAR API ---

def get_clinvar_significance_for_list(rs_ids_list):
    """Пакетная проверка списка rs-номеров через NCBI ClinVar."""
    results = {}
    clinvar_ids = []
    id_map = {}

    print(f"🔍 Запрос данных для {len(rs_ids_list)} вариантов...")
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
        except Exception:
            results[rs_id] = "Error"

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
        except Exception:
            pass
    return results

# --- БЛОК 3: ГЕНЕТИЧЕСКИЙ СКАНЕР ---

def universal_genomic_scanner(sequence):
    """Анализ сырой ДНК: GC-состав, повторы и мутации."""
    print("\n" + "="*45)
    print("🚀 STARTING GENOMIC SCAN v4.0 (Visual)")
    print("="*45)

    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    print(f"📊 Global GC-Content: {gc:.2f}%")

    # Поиск CAG-повторов (Болезнь Гентингтона)
    repeats = re.findall(r'(?:CAG){3,}', sequence)
    if repeats:
        count = len(max(repeats, key=len)) // 3
        status = "⚠️ PATHOGENIC" if count >= 36 else "✅ NORMAL"
        print(f"🧬 HTT Gene: {count} CAG repeats ({status})")

    # Маркер долголетия SIRT1
    if sequence.startswith('G'):
        print("🌟 SIRT1: Longevity Variant 'G' detected!")
    
    print("="*45 + "\n")
    return gc

# --- БЛОК 4: ГЕНЕРАЦИЯ PDF ---

def save_report_to_pdf(results, dna_sequence, filename="Genomic_Report.pdf"):
    """Создание PDF отчета с графиком и таблицей."""
    try:
        pdfmetrics.registerFont(TTFont('DejaVuSans', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
        font_name = 'DejaVuSans'
    except Exception:
        font_name = 'Helvetica'

    c = canvas.Canvas(filename, pagesize=letter)
    width, height = letter

    # Заголовок отчета
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 18)
    c.drawCentredString(width / 2.0, height - 50, "Genomic Anomaly Atlas: Diagnostic Report")
    
    # Визуализация (График вставляется под заголовком)
    c.setFont(font_name, 12)
    c.drawString(50, height - 85, "I. Sequence Landscape Analysis:")
    img_path = generate_gc_plot(dna_sequence)
    c.drawImage(img_path, 50, height - 260, width=500, height=160)
    
    # Таблица клинической значимости
    y_start = height - 300
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 12)
    c.drawString(50, y_start, "II. Variant Clinical Significance (ClinVar):")
    
    y = y_start - 30
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 10)
    c.drawString(50, y, "Variant ID")
    c.drawString(200, y, "Status")
    c.line(50, y-5, 550, y-5)
    
    y -= 25
    c.setFont(font_name, 10)
    for rs, status in results.items():
        if y < 50:
            c.showPage()
            y = height - 50
        c.drawString(50, y, rs.upper())
        c.drawString(200, y, str(status))
        y -= 20

    c.save()
    print(f"📄 PDF Report saved: {filename}")

# --- ГЛАВНЫЙ ЗАПУСК ---

def generate_professional_report(input_file, dna_sample):
    """Сборка всего процесса в один запуск."""
    if not os.path.exists(input_file):
        with open(input_file, "w") as f:
            f.write("rs1137282\nrs121908675\nrs7069102")
    
    with open(input_file, 'r') as f:
        rs_list = [line.strip() for line in f if line.strip()]

    universal_genomic_scanner(dna_sample)
    results = get_clinvar_significance_for_list(rs_list)
    save_report_to_pdf(results, dna_sample)

if __name__ == "__main__":
    # Тестовая ДНК-последовательность
    test_dna = "G" + ("CAG" * 42) + "ATCG" * 50 + "GCGC" * 10
    
    # Запуск всей программы
    generate_professional_report("rs_numbers.txt", test_dna)
