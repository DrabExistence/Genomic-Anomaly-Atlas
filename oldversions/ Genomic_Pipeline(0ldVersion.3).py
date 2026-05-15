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
    plt.title("Геномный Ландшафт: Распределение GC-состава", fontsize=12) # Translated
    plt.xlabel("Позиция в Последовательности (bp)", fontsize=10) # Translated
    plt.ylabel("GC %", fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_img, dpi=100)
    plt.close()
    return output_img

# --- БЛОК 2: РАБОТА С CLINVAR API ---

def get_clinvar_significance_for_list(rs_ids_list):
    """Запрос данных из NCBI ClinVar."""
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
                id_map[cid] = rs_id # Store mapping
            else:
                results[rs_id] = "Не найдено" # Explicitly mark as Not Found, translated
        except Exception as e:
            results[rs_id] = f"Ошибка API (поиск ClinVar): {e}" # More specific error, translated

    # After initial esearch, any rs_id not in results yet is a candidate for batch esummary
    # rs_ids_for_esummary = {id_map[cid] for cid in clinvar_ids} # Not directly used now, but useful for understanding

    if clinvar_ids:
        try:
            handle = Entrez.esummary(db="clinvar", id=",".join(clinvar_ids), validate=False)
            summaries = Entrez.read(handle, validate=False)
            handle.close()
            processed_rs_ids = set() # Keep track of rs_ids successfully processed by esummary

            for doc in summaries['DocumentSummarySet']['DocumentSummary']:
                uid = str(doc.get('uid') or getattr(doc, 'attributes', {}).get('uid', ''))
                orig_rs = id_map.get(uid)
                if orig_rs:
                    results[orig_rs] = doc['germline_classification']['description']
                    processed_rs_ids.add(orig_rs)

            # Identify rs_ids that found a ClinVar ID but failed to get summary
            for cid in clinvar_ids:
                orig_rs = id_map.get(cid)
                if orig_rs not in processed_rs_ids and orig_rs not in results: # If not processed and not already marked as Not Found from esearch
                    results[orig_rs] = "Ошибка API (нет описания ClinVar)" # Explicit error for failed summary, translated

        except Exception as e:
            # If batch esummary itself fails, mark all pending rs_ids as API Error
            for cid in clinvar_ids:
                orig_rs = id_map.get(cid)
                if orig_rs not in results: # Only if not already marked
                    results[orig_rs] = f"Ошибка API (сбой пакетного запроса ClinVar): {e}" # Translated

    # Ensure all original rs_ids in rs_ids_list are accounted for in results
    for rs_id in rs_ids_list:
        if rs_id not in results:
            results[rs_id] = "Ошибка API (неизвестная ошибка ClinVar)" # Fallback for any unhandled case, translated

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
    # Default fonts
    font_regular_current = 'Helvetica'
    font_bold_current = 'Helvetica-Bold'

    try:
        dejavu_sans_path = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
        dejavu_sans_bold_path = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'

        if os.path.exists(dejavu_sans_path) and os.path.exists(dejavu_sans_bold_path):
            # Register with distinct names to avoid conflicts
            pdfmetrics.registerFont(TTFont('DejaVuSansRegular', dejavu_sans_path))
            pdfmetrics.registerFont(TTFont('DejaVuSansBold', dejavu_sans_bold_path))
            font_regular_current = 'DejaVuSansRegular'
            font_bold_current = 'DejaVuSansBold'
            print(f"✅ Шрифты DejaVuSans успешно зарегистрированы.") # Translated
        else:
            print("⚠️ Внимание: Файлы шрифтов DejaVuSans не найдены. Используется Helvetica.") # Translated
    except Exception as e:
        print(f"❌ Ошибка при регистрации шрифтов DejaVuSans: {e}. Используется Helvetica.") # Translated

    c = canvas.Canvas(filename, pagesize=letter)
    width, height = letter

    c.setFont(font_bold_current, 18)
    c.drawCentredString(width / 2.0, height - 50, "Геномный Атлас Аномалий: Диагностический Отчет") # Translated

    c.setFont(font_regular_current, 12)
    c.drawString(50, height - 85, "I. Анализ Ландшафта Последовательности:") # Translated
    img_path = generate_gc_plot(dna_sequence)
    c.drawImage(img_path, 50, height - 260, width=500, height=160)

    curr_y = height - 280
    c.setFont(font_bold_current, 12)
    c.drawString(50, curr_y, "II. Основные Метрики Сканирования:") # Translated
    c.setFont(font_regular_current, 10)
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
    c.setFont(font_bold_current, 12)
    c.drawString(50, curr_y, "III. Анализ Сходства Экзонов:") # Translated
    c.setFont(font_regular_current, 10)
    for name, sim in exon_results.items():
        curr_y -= 20
        c.drawString(60, curr_y, f"• {name}: {sim}% сходства с образцом") # Translated

    curr_y -= 35
    c.setFont(font_bold_current, 12)
    c.drawString(50, curr_y, "IV. Клиническая Значимость (ClinVar):") # Translated
    curr_y -= 25
    c.setFont(font_regular_current, 9)

    if "INFO_MESSAGE" in results: # Added INFO_MESSAGE handling
        message = results["INFO_MESSAGE"]
        try:
            from reportlab.lib.enums import TA_LEFT
            from reportlab.platypus import Paragraph
            from reportlab.lib.styles import getSampleStyleSheet

            styles = getSampleStyleSheet()
            style = styles['Normal']
            style.fontName = font_regular_current
            style.fontSize = 9 # Adjusted for ClinVar table
            style.alignment = TA_LEFT

            p = Paragraph(message, style)
            w, h = p.wrapOn(c, width - 100, height)
            p.drawOn(c, 50, curr_y - h)
            curr_y -= h + 12

        except ImportError:
            # Fallback if reportlab.platypus is not available
            wrapped_lines = wrap_text(message, font_regular_current, 9, width - 100, c)
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
                c.setFont(font_regular_current, 9)
                c.drawString(60, curr_y, "ID Варианта") # Translated
                c.drawString(180, curr_y, "Статус") # Translated
                curr_y -= 10
                c.line(50, curr_y, width - 50, curr_y)
                curr_y -= 15

            display_status = "Запись ClinVar не найдена" if status == "Не найдено" else str(status) # Translated
            c.drawString(60, curr_y, rs.upper())
            c.drawString(180, curr_y, display_status)
            curr_y -= 15

    c.save()
    print(f"✅ Отчет '{filename}' сохранен локально.") # Translated
    if os.path.exists("gc_plot.png"): os.remove("gc_plot.png")
    return filename # Return the filename

# --- НОВЫЙ БЛОК: ПАРСЕР VCF ---

def parse_vcf(vcf_path):
    """Извлекает RS-номера из VCF файла."""
    rs_numbers = []
    if not os.path.exists(vcf_path):
        print(f"❌ Файл {vcf_path} не найден.")
        return []
    print(f"📂 Парсинг VCF файла: {vcf_path}...") # Translated

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

        print(f"✅ Найдено {len(rs_numbers)} вариантов в VCF.") # Translated
        return rs_numbers
    except Exception as e:
        print(f"❌ Ошибка при чтении VCF: {e}") # Translated
        return []

# --- Вспомогательные функции для CLI ---
def _load_file_content(file_path, default_content, warning_message):
    """Загружает содержимое файла или возвращает значение по умолчанию."""
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            return f.read().strip()
    else:
        print(warning_message)
        return default_content

def _get_rs_list_and_clinvar_info(args):
    """Извлекает список RS-номеров и формирует начальное сообщение ClinVar.
    Возвращает (rs_list, clinvar_data).
    """
    rs_list = []
    clinvar_data = {}

    if args.vcf: # If VCF file is specified
        rs_list = parse_vcf(args.vcf)
        if not rs_list:
            clinvar_data["INFO_MESSAGE"] = "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. VCF файл обработан, но RS-номера не найдены или произошла ошибка."
    elif args.input: # If RS numbers file is specified
        if os.path.exists(args.input):
            with open(args.input, 'r') as f:
                rs_list = [l.strip() for l in f if l.strip()]
            if not rs_list:
                clinvar_data["INFO_MESSAGE"] = "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. Файл с RS-номерами пуст."
        else:
            print(f"❌ Файл {args.input} не найден. Невозможно получить RS-номера.")
            clinvar_data["INFO_MESSAGE"] = f"Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. Файл '{args.input}' не найден."
    else: # Neither --vcf nor --input was provided, use default RS numbers
        print("⚠️ Ни VCF файл, ни файл с RS-номерами не предоставлены. Использую RS-номера по умолчанию.")
        rs_list = ["rs1137282", "rs121908675", "rs7069102"]
        clinvar_data["INFO_MESSAGE"] = "Раздел 'Клиническая Значимость (ClinVar)' не может быть сгенерирован. Ни VCF файл, ни файл с RS-номерами не предоставлены; использованы RS-номера по умолчанию."

    return rs_list, clinvar_data

# --- БЛОК 5: CLI ИНТЕРФЕЙС ---

def main():
    parser = argparse.ArgumentParser(
        description="🧬 Геномный Атлас Аномалий: Диагностический Пайплайн v4.5", # Version updated
        epilog="Пример: python Genomic_Pipeline.py --vcf test_data.vcf -d dna.txt -o Report.pdf\n        ИЛИ: python Genomic_Pipeline.py -i rs_numbers.txt -d dna.txt -o Report.pdf"
    )

    # Создаем взаимоисключающую группу для источников RS-номеров (файл или VCF)
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "-i", "--input",
        help="Путь к файлу с RS-номерами (например, rs_numbers.txt). Не может использоваться с --vcf."
    )
    input_group.add_argument(
        "-v", "--vcf",
        help="Путь к VCF файлу для анализа. Не может использоваться с --input."
    )

    parser.add_argument("-d", "--dna", default="sample_dna.txt", help="Путь к файлу ДНК (по умолчанию: sample_dna.txt).") # Translated
    parser.add_argument("-o", "--output", default="Genomic_Report.pdf", help="Имя выходного PDF файла (по умолчанию: Genomic_Report.pdf).") # Translated

    args = parser.parse_args() # Используем parse_args() для строгой проверки аргументов

    # 1. Получение ДНК
    dna_sample = _load_file_content(
        args.dna,
        "G" + ("CAG" * 42) + "ATCG" * 100,
        f"⚠️ Файл {args.dna} не найден. Использую тестовую ДНК."
    )

    # 2. Получение списка RS и начального сообщения ClinVar
    rs_list, clinvar_data = _get_rs_list_and_clinvar_info(args)

    # 3. Запрос данных ClinVar, если список RS не пуст и нет INFO_MESSAGE
    if not clinvar_data and rs_list:
        clinvar_data = get_clinvar_significance_for_list(rs_list)

    print("🚀 Запуск анализа...") # Translated
    scan_results = universal_genomic_scanner(dna_sample)

    exon_8_ref = "ATTATAGTAGAGAGA"
    sim_8 = calculate_similarity(dna_sample[:len(exon_8_ref)], exon_8_ref)

    exon_9_ref = "ATTAGAGTAGAGACA"
    sim_9 = calculate_similarity(dna_sample[:len(exon_9_ref)], exon_9_ref)

    exon_comparison = {"Exon_8_IIIb_Reference": sim_8, "Exon_9_IIIc_Reference": sim_9}

    # save_report_to_pdf now returns the filename
    final_report_name = save_report_to_pdf(clinvar_data, dna_sample, scan_results, exon_comparison, args.output)

    # Conditional download logic for Colab
    try:
        from google.colab import files
        if os.path.exists(final_report_name):
            files.download(final_report_name)
            print(f"✅ Отчет '{final_report_name}' успешно загружен.")
        else:
            print(f"⚠️ Файл отчета '{final_report_name}' не найден для загрузки.")
    except ImportError:
        print("ℹ️ google.colab не найдена. Файл сохранен локально, загрузка не требуется.")
    except Exception as e:
        print(f"❌ Ошибка при попытке загрузки файла: {e}")

    print(f"✅ Готово! Отчет сохранен в: {final_report_name}") # Translated
    return final_report_name # Return the filename from main

if __name__ == "__main__":
    main()
