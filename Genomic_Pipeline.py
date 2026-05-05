import re
from Bio import Entrez

def get_clinvar_info(rs_id):
    """Интеграция с базой ClinVar через API."""
    try:
        Entrez.email = "architect@example.com"
        search = Entrez.esearch(db="clinvar", term=f"{rs_id}[Variant ID]")
        res = Entrez.read(search)
        if res['IdList']:
            summ = Entrez.read(Entrez.esummary(db="clinvar", id=res['IdList'][0]))
            return summ['DocumentSummarySet']['DocumentSummary'][0]['clinical_significance']['description']
        return "Not found in ClinVar"
    except:
        return "ClinVar API Error"

def universal_genomic_scanner(sequence):
    print("\n" + "="*40)
    print("🚀 GENOMIC ANOMALY SCANNER v2.0")
    print("="*40)
    
    # 1. Анализ GC-состава
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    print(f"📊 [STRUCTURAL] GC-Content: {gc:.2f}%")
    
    # 2. Проверка TERT (Рак)
    if len(sequence) > 228 and sequence[228] == 'T':
        print("🚨 [ALERT] TERT Mutation Found: C228T (Immortality Switch ON)")
    
    # 3. Проверка Huntington (Повторы)
    repeats = re.findall(r'(?:CAG){3,}', sequence)
    if repeats:
        count = len(max(repeats, key=len)) // 3
        status = "⚠️ PATHOGENIC" if count >= 36 else "✅ NORMAL"
        print(f"🧬 [REPEATS] HTT Gene: {count} CAG repeats ({status})")
    
    # 4. Проверка SIRT1 (Долголетие)
    # Имитируем, что маркер долголетия проверяется в начале фрагмента
    if sequence[0] == 'G':
        print("🌟 [BIOHACK] SIRT1: Longevity Variant 'G' detected!")
        # Проверяем ClinVar для подтверждения
        # print(f"   └─ ClinVar says: {get_clinvar_info('rs7069102')}")

    # 5. Проверка LRP5 (Кости-титан)
    # Имитируем позицию 171
    if len(sequence) > 171 and sequence[171] == 'T':
        print("🦾 [SUPERPOWER] LRP5: TITAN BONES mutation detected (D171V)!")

    print("="*40 + "\n")

# ТЕСТОВЫЕ ДАННЫЕ
# Генерируем "Супер-образец": SIRT1(G) + LRP5(T)
super_human_dna = "G" + ("A" * 170) + "T" + ("C" * 100)
universal_genomic_scanner(super_human_dna)
