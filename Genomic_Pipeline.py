import re  # Стандартная библиотека для поиска текста
from Bio import Entrez  # Библиотека для работы с биологическими БД
from clinvar_api import get_clinvar_significance_for_list, load_rs_from_file  # новый модуль связи с ClinVar
from fpdf import FPDF

def universal_genomic_scanner(sequence):
    print("\n" + "="*45)
    print("🚀 STARTING MULTI-LEVEL GENOMIC SCAN v3.0")
    print("="*45)
    
    # 1. Расчет GC-состава
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    print(f"📊 [STRUCTURAL] GC-Content: {gc:.2f}%")
    
    # 2. ПРОВЕРКА ПОВТОРОВ (HTT)
    repeats = re.findall(r'(?:CAG){3,}', sequence)
    if repeats:
        count = len(max(repeats, key=len)) // 3
        status = "⚠️ PATHOGENIC" if count >= 36 else "✅ NORMAL"
        print(f"🧬 [REPEATS] HTT Gene: {count} CAG repeats ({status})")
    else:
        print(f"🧬 [REPEATS] HTT Gene: No expansions detected (Normal)")

    # 3. ПРОВЕРКА SIRT1 (Долголетие)
    if sequence[0] == 'G':
        print("🌟 [BIOHACK] SIRT1: Longevity Variant 'G' detected!")
        print("🔍 Connecting to ClinVar to verify...")
        clinical_status = get_clinvar_significance('rs7069102') 
        print(f"   └─ Official Scientific Status: {clinical_status}")

    # 4. ПРОВЕРКА LRP5 (Кости-титан)
    if len(sequence) > 171 and sequence[171] == 'T':
        print("🦾 [SUPERPOWER] LRP5: High Bone Density variant found!")
        print("🔍 Verifying Clinical Significance...")
        clinical_status = get_clinvar_significance('rs121908675')
        print(f"   └─ ClinVar Result: {clinical_status}")

    # 5. ПРОВЕРКА TERT (Кнопка рака)
    if len(sequence) > 228 and sequence[228] == 'T':
        print("🚨 [CANCER RISK] TERT Promoter Mutation Found (C228T)")

    print("="*45 + "\n")

# ТЕСТОВЫЕ ДАННЫЕ (Полная проверка всех систем)
# Добавили CAG-повторы, чтобы пункт 2 тоже сработал
super_human_dna = "G" + ("CAG" * 40) + ("A" * 48) + "T" + ("A" * 56) + "T" + ("C" * 20)
universal_genomic_scanner(super_human_dna)

def save_report_to_pdf(results, filename="Genomic_Report.pdf"):
    """Создает PDF-файл на основе словаря с результатами."""
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(200, 10, txt="Genomic Anomaly Atlas: Diagnostic Report", ln=True, align='C')
    pdf.ln(10)
    
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(50, 10, "Variant (ID)", border=1)
    pdf.cell(140, 10, "Clinical Significance", border=1)
    pdf.ln()
    
    pdf.set_font("Arial", size=10)
    for rs, status in results.items():
        pdf.cell(50, 10, txt=rs.upper(), border=1)
        pdf.cell(140, 10, txt=str(status), border=1)
        pdf.ln()
    
    pdf.output(filename)
    print(f"📄 PDF Report saved as: {filename}")

def generate_professional_report(input_file):
    """
    Массовая проверка мутаций из файла и автоматическое создание PDF.
    """
    print("\n" + "═"*60)
    print("📋 GENOMIC ANOMALY ATLAS: BATCH PROCESSING")
    print("═"*60)
    
    # 1. Загружаем список rs-номеров из файла
    rs_list = load_rs_from_file(input_file)
    if not rs_list:
        print("❌ Error: No variants found in the input file.")
        return
    
    # 2. Получаем данные из ClinVar (пакетный запрос)
    print(f"🔍 Fetching data for {len(rs_list)} variants...")
    results = get_clinvar_significance_for_list(rs_list)
    
    # 3. Вывод в консоль для быстрой проверки
    for rs, status in results.items():
        print(f" -> {rs.upper()}: {status}")
    
    # 4. ВЫЗОВ PDF ГЕНЕРАТОРА (та самая 'основная логика')
    save_report_to_pdf(results)
    
    print("═"*60 + "\n")

# --- ГЛАВНЫЙ ЗАПУСК ---
if __name__ == "__main__":
    # Тест сканера сырой ДНК
    test_dna = "G" + ("CAG" * 40) + ("A" * 48) + "T" + ("A" * 56) + "T" + ("C" * 20)
    universal_genomic_scanner(test_dna)
    
    # Тест пакетной обработки и генерации PDF
    # Убедись, что файл rs_numbers.txt лежит в той же папке в Colab!
    generate_professional_report("rs_numbers.txt")
