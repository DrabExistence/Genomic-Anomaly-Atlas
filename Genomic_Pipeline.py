import re
import os
from Bio import Entrez
from tqdm.auto import tqdm
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# --- НАСТРОЙКИ NCBI ---
Entrez.email = "zzzender20@gmail.com" 

# --- БЛОК 1: РАБОТА С CLINVAR API ---

def get_clinvar_significance(rs_id):
    """Получает клиническую значимость для одного rs-номера."""
    try:
        search_query = f"{rs_id}[Variant ID]"
        handle = Entrez.esearch(db="clinvar", term=search_query)
        search_results = Entrez.read(handle)
        handle.close()
        
        if not search_results['IdList']:
            return "No ClinVar record found."

        clinvar_id = search_results['IdList'][0]
        handle = Entrez.esummary(db="clinvar", id=clinvar_id)
        summary = Entrez.read(handle, validate=False)
        handle.close()

        return summary['DocumentSummarySet']['DocumentSummary'][0]['germline_classification']['description']
    except Exception:
        return "Data unavailable."

def get_clinvar_significance_for_list(rs_ids_list):
    """Пакетная проверка списка rs-номеров."""
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
        except:
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
        except:
            pass
    return results

def load_rs_from_file(filepath):
    """Загрузка rs-номеров из файла."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# --- БЛОК 2: ГЕНЕТИЧЕСКИЙ СКАНЕР ---

def universal_genomic_scanner(sequence):
    """Анализ сырой последовательности ДНК."""
    print("\n" + "="*45)
    print("🚀 STARTING GENOMIC SCAN v3.0")
    print("="*45)

    # GC-состав
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    print(f"📊 GC-Content: {gc:.2f}%")

    # CAG Повторы (Гентингтон)
    repeats = re.findall(r'(?:CAG){3,}', sequence)
    if repeats:
        count = len(max(repeats, key=len)) // 3
        status = "⚠️ PATHOGENIC" if count >= 36 else "✅ NORMAL"
        print(f"🧬 HTT Gene: {count} CAG repeats ({status})")

    # SIRT1 (Долголетие)
    if sequence.startswith('G'):
        print("🌟 SIRT1: Longevity Variant 'G' detected!")
    
    print("="*45 + "\n")

# --- БЛОК 3: ГЕНЕРАЦИЯ ОТЧЕТА (PDF) ---

def save_report_to_pdf(results, filename="Genomic_Report.pdf"):
    """Создание PDF отчета с поддержкой кириллицы."""
    try:
        # Пытаемся подключить шрифт DejaVu (стандарт для Linux/Colab)
        pdfmetrics.registerFont(TTFont('DejaVuSans', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
        font_name = 'DejaVuSans'
    except:
        # Если шрифта нет (например, на Windows), используем стандартный Helvetica
        font_name = 'Helvetica'

    c = canvas.Canvas(filename, pagesize=letter)
    width, height = letter

    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 16)
    c.drawCentredString(width / 2.0, height - 50, "Genomic Anomaly Atlas: Diagnostic Report")

    y = height - 100
    c.setFont(font_name, 12)
    c.drawString(50, y, "Variant (ID)")
    c.drawString(200, y, "Clinical Significance")
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

def generate_professional_report(input_file):
    """Основной процесс: Файл -> API -> PDF."""
    rs_list = load_rs_from_file(input_file)
    if not rs_list:
        # Создаем тестовый файл, если его нет
        with open(input_file, "w") as f:
            f.write("rs1137282\nrs121908675\nrs7069102")
        rs_list = load_rs_from_file(input_file)

    results = get_clinvar_significance_for_list(rs_list)
    save_report_to_pdf(results)

# --- ЗАПУСК ---
if __name__ == "__main__":
    # 1. Анализ тестовой ДНК
    dna_sample = "G" + ("CAG" * 40) + "ATCG" * 20
    universal_genomic_scanner(dna_sample)
    
    # 2. Пакетный отчет
    generate_professional_report("rs_numbers.txt")
