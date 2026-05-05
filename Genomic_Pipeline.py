import re  # Стандартная библиотека для поиска текста
from Bio import Entrez  # Библиотека для работы с биологическими БД
from clinvar_api import get_clinvar_significance  # новый модуль связи с ClinVar

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
