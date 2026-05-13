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
from google.colab import files # Added for files.download in save_report_to_pdf
import io

# --- НАСТРОЙКИ NCBI ---
Entrez.email = "genomic.atlas.project@example.com"

# --- БЛОК 1: МАТЕМАТИКА И ВИЗУАЛИЗАЦИЯ ---

def calculate_similarity(seq1, seq2):
    min_len = min(len(seq1), len(seq2))
    if max(len(seq1), len(seq2)) == 0: return 0
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)

    # Сходство рассчитывается на основе максимальной длины, чтобы учесть пропуски/вставки
    similarity = (matches / max(len(seq1), len(seq2))) * 100
    return round(similarity, 2)

def generate_gc_plot(sequence, window_size=20, output_img="gc_plot.png"):
    values = []
    for i in range(len(sequence) - window_size + 1):
        subseq = sequence[i:i + window_size]
        gc = (subseq.count('G') + subseq.count('C')) / len(subseq) * 100
        values.append(gc)

    plt.figure(figsize=(8, 3))
    plt.plot(values, color='#2c3e50', linewidth=1.5, label='GC %')
    plt.fill_between(range(len(values)), values, color='#3498db', alpha=0.3)
    plt.axhline(y=50, color='r', linestyle='--', alpha=0.3)
    plt.title("Геномный Ландшафт: Радар GC-Состава", fontsize=12) # Translated
    plt.xlabel("Позиция в Последовательности (bp)", fontsize=10) # Translated
    plt.ylabel("GC %", fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_img, dpi=100)
    plt.close()
    return output_img

# --- БЛОК 2: РАБОТА С CLINVAR API ---

def get_clinvar_significance_for_list(rs_ids_list):
    results = {}
    clinvar_ids = []
    id_map = {}

    if not rs_ids_list:
        return results

    print(f"🔍 Запрос данных ClinVar для {len(rs_ids_list)} вариантов...")
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
            results[rs_id] = "Connection Error"

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
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

    repeats = re.findall(r'(?:CAG){3,}', sequence)
    htt_count = len(max(repeats, key=len)) // 3 if repeats else 0
    htt_status = "PATHOGENIC" if htt_count >= 36 else "NORMAL"

    sirt1_detected = sequence.startswith('G')

    return {
        "gc_content": gc,
        "htt_count": htt_count,
        "htt_status": htt_status,
        "sirt1": sirt1_detected
    }

def wrap_text(text, font_name, font_size, max_width, canvas_obj):
    paragraphs = text.split('\n')
    all_wrapped_lines = []

    for paragraph in paragraphs:
        if not paragraph.strip():
            all_wrapped_lines.append('')
            continue

        words = paragraph.split(' ')
        current_line = []
        for word in words:
            test_line = ' '.join(current_line + [word])
            if canvas_obj.stringWidth(test_line, font_name, font_size) < max_width:
                current_line.append(word)
            else:
                if current_line:
                    all_wrapped_lines.append(' '.join(current_line))
                current_line = [word]
        if current_line:
            all_wrapped_lines.append(' '.join(current_line))
    return all_wrapped_lines

# --- БЛОК 4: ГЕНЕРАЦИЯ ОТЧЕТА ---

def save_report_to_pdf(results, dna_sequence, scan_data, exon_results, filename):
    try:
        pdfmetrics.registerFont(TTFont('DejaVuSans', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
        pdfmetrics.registerFont(TTFont('DejaVuSans-Bold', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf')) # Added bold font
        font_name = 'DejaVuSans'
    except:
        font_name = 'Helvetica'

    c = canvas.Canvas(filename, pagesize=letter)
    width, height = letter

    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 18)
    c.drawCentredString(width / 2.0, height - 50, "Геномный Атлас Аномалий: Диагностический Отчет") # Translated

    c.setFont(font_name, 12)
    c.drawString(50, height - 85, "I. Анализ Ландшафта Последовательности:") # Translated
    img_path = generate_gc_plot(dna_sequence)
    c.drawImage(img_path, 50, height - 260, width=500, height=160)

    curr_y = height - 280
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 12)
    c.drawString(50, curr_y, "II. Основные Метрики Сканирования:") # Translated
    c.setFont(font_name, 10)
    curr_y -= 20
    c.drawString(60, curr_y, f"• Глобальное содержание GC: {scan_data['gc_content']:.2f}%") # Translated
    curr_y -= 15

    htt_status_ru = "" # Prepare for translation
    if scan_data['htt_status'] == "PATHOGENIC":
        htt_status_ru = "ПАТОГЕННЫЙ"
    elif scan_data['htt_status'] == "NORMAL":
        htt_status_ru = "НОРМАЛЬНЫЙ"
    c.drawString(60, curr_y, f"• Повторы HTT: {scan_data['htt_count']} ({htt_status_ru})") # Translated
    curr_y -= 15

    sirt1_status_ru = "ДА" if scan_data['sirt1'] else "НЕТ" # Prepare for translation
    c.drawString(60, curr_y, f"• Маркер долголетия SIRT1: {sirt1_status_ru}") # Translated

    curr_y -= 30
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 12)
    c.drawString(50, curr_y, "III. Анализ Сходства Экзонов:") # Translated
    c.setFont(font_name, 10)
    for name, sim in exon_results.items():
        curr_y -= 20
        c.drawString(60, curr_y, f"• {name}: {sim}% сходства с образцом") # Translated

    curr_y -= 35
    c.setFont(f"{font_name}-Bold" if font_name != 'Helvetica' else 'Helvetica-Bold', 12)
    c.drawString(50, curr_y, "IV. Клиническая Значимость (ClinVar):") # Translated
    curr_y -= 25
    c.setFont(font_name, 9)

    if "INFO_MESSAGE" in results: # Added INFO_MESSAGE handling
        message = results["INFO_MESSAGE"]
        try:
            from reportlab.lib.enums import TA_LEFT
            from reportlab.platypus import Paragraph
            from reportlab.lib.styles import getSampleStyleSheet

            styles = getSampleStyleSheet()
            style = styles['Normal']
            style.fontName = font_name
            style.fontSize = 9 # Adjusted for ClinVar table
            style.alignment = TA_LEFT

            p = Paragraph(message, style)
            w, h = p.wrapOn(c, width - 100, height)
            p.drawOn(c, 50, curr_y - h)
            curr_y -= h + 12

        except ImportError:
            wrapped_lines = wrap_text(message, font_name, 9, width - 100, c)
            for line in wrapped_lines:
                c.drawString(50, curr_y, line)
                curr_y -= 12
    else:
        # Headers for ClinVar table
        c.drawString(60, curr_y, "ID Варианта") # Translated
        c.drawString(180, curr_y, "Статус") # Translated
        curr_y -= 10
        c.line(50, curr_y, width - 50, curr_y) # Line under headers
        curr_y -= 15

        for rs, status in results.items():
            if curr_y < 50:
                c.showPage()
                curr_y = height - 50
                c.setFont(font_name, 9)
                c.drawString(60, curr_y, "ID Варианта") # Translated
                c.drawString(180, curr_y, "Статус") # Translated
                curr_y -= 10
                c.line(50, curr_y, width - 50, curr_y)
                curr_y -= 15

            display_status = "Запись ClinVar не найдена" if status == "Not Found" else str(status) # Translated
            c.drawString(60, curr_y, rs.upper())
            c.drawString(180, curr_y, display_status)
            curr_y -= 15

    c.save()
    # files.download(filename) # REMOVED THIS LINE
    if os.path.exists("gc_plot.png"): os.remove("gc_plot.png")

# --- НОВЫЙ БЛОК: ПАРСЕР VCF ---

def parse_vcf(vcf_path):
    """Извлекает RS-номера из VCF файла."""
    rs_numbers = []
    print(f"📂 Парсинг VCF файла: {vcf_path}...")

    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                # Пропускаем заголовки (начинаются с #)
                if line.startswith('#'):
                    continue

                # VCF колонки: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
                columns = line.split('\t')
                if len(columns) > 2:
                    variant_id = columns[2] # Колонка ID обычно содержит rs...
                    if variant_id.startswith('rs'):
                        rs_numbers.append(variant_id)

        print(f"✅ Найдено {len(rs_numbers)} вариантов в VCF.")
        return rs_numbers
    except Exception as e:
        print(f"❌ Ошибка при чтении VCF: {e}")
        return []

# --- БЛОК 5: CLI ИНТЕРФЕЙС ---

def main():
    parser = argparse.ArgumentParser(
        description="🧬 Геномный Атлас Аномалий: Диагностический Пайплайн v4.0", # Translated
        epilog="Пример: python Genomic_Pipeline.py -i rs_numbers.txt -d dna.txt -o Report.pdf" # Epilog is fine as example
    )
    parser.add_argument("-i", "--input", default="rs_numbers.txt", help="Файл с RS-номерами") # Translated
    parser.add_argument("-v", "--vcf", help="Путь к VCF файлу для анализа") # НОВЫЙ АРГУМЕНТ
    parser.add_argument("-d", "--dna", default="sample_dna.txt", help="Файл с ДНК последовательностью") # Translated
    parser.add_argument("-o", "--output", default="Genomic_Report.pdf", help="Имя выходного PDF") # Translated

    args, unknown = parser.parse_known_args()

    dna_sample = ""
    if not os.path.exists(args.dna):
        print(f"⚠️ Файл {args.dna} не найден. Использую тестовую ДНК.") # Translated
        dna_sample = "G" + ("CAG" * 42) + "ATCG" * 100
    else:
        with open(args.dna, 'r') as f: dna_sample = f.read().strip()

    rs_list = []
    clinvar_data = {}
    if args.vcf:
        rs_list = parse_vcf(args.vcf)
        if not rs_list:
            clinvar_data = {"INFO_MESSAGE": "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. VCF файл обработан, но RS-номера не найдены или произошла ошибка."}
    elif not os.path.exists(args.input):
        print(f"⚠️ Файл {args.input} не найден.") # Translated
        clinvar_data = {"INFO_MESSAGE": "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. Файл 'rs_numbers.txt' не найден или оказался пустым. Для получения результатов, пожалуйста, создайте файл 'rs_numbers.txt' в корневой директории с одним rsID гена на каждой новой строке."}
    else:
        with open(args.input, 'r') as f:
            rs_list = [l.strip() for l in f if l.strip()]
        if not rs_list:
            print(f"⚠️ Файл {args.input} пуст.") # Translated
            clinvar_data = {"INFO_MESSAGE": "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. Файл 'rs_numbers.txt' не найден или оказался пустым. Для получения результатов, пожалуйста, создайте файл 'rs_numbers.txt' в корневой директории с одним rsID гена на каждой новой строке."}

    if not clinvar_data and rs_list:
        clinvar_data = get_clinvar_significance_for_list(rs_list)
    elif not clinvar_data and not rs_list and not args.vcf and not os.path.exists(args.input):
        # Default message if no input source provided and no default rs_numbers.txt
        clinvar_data = {"INFO_MESSAGE": "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. Не предоставлен ни VCF файл, ни файл с RS-номерами."}

    print("🚀 Запуск анализа...") # Translated
    scan_results = universal_genomic_scanner(dna_sample)

    exon_8_ref = "ATTATAGTAGAGAGA"
    sim_8 = calculate_similarity(dna_sample[:len(exon_8_ref)], exon_8_ref)

    exon_9_ref = "ATTAGAGTAGAGACA"
    sim_9 = calculate_similarity(dna_sample[:len(exon_9_ref)], exon_9_ref)

    exon_comparison = {"Exon_8_IIIb_Reference": sim_8, "Exon_9_IIIc_Reference": sim_9}

    save_report_to_pdf(clinvar_data, dna_sample, scan_results, exon_comparison, args.output)
    print(f"✅ Готово! Отчет сохранен в: {args.output}") # Translated

if __name__ == "__main__":
    main()
