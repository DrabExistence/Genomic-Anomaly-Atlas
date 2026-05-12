# AI Pipelines and Automation (AI-Driven Pipeline)
Now we'll combine all these skills into a single automation tool to create a "pipeline" that takes raw data and produces a finished report.

# What's the idea?
Instead of running five different scripts, we'll create a Master Script. It will:

Check the GC composition (to determine if it's a promoter).

Look for mutations (like in TERT).

Check for repeats (like in HTT).

In this file we will describe the logic of the "Universal Analyzer".

Python
```python
import re

def universal_genomic_scanner(sequence):
    print("--- Starting Genomic Analysis Pipeline ---")
    
    # 1. GC-Content (Работает отлично!)
    gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    print(f"[STEP 1] GC-Content: {gc:.2f}%")
    
    # 2. Поиск мутации TERT (Исправлено)
    # Проверяем 228-й нуклеотид (помним, что в Python счет с 0)
    if len(sequence) > 228:
        if sequence[228] == 'T':
            print("[STEP 2] ALERT: TERT Promoter Mutation (C228T) detected! Immortality switch ON.")
    
    # 3. Поиск повторов CAG (Работает отлично!)
    repeats = re.findall(r'(?:CAG){3,}', sequence)
    if repeats:
        count = len(max(repeats, key=len)) // 3
        if count >= 36:
            print(f"[STEP 3] PATHOGENIC REPEAT: Found {count} CAG repeats (Huntington Risk)")
        else:
            print(f"[STEP 3] Normal repeats: {count} CAG")
        
    print("--- Analysis Complete ---")

# Твой пайплайн готов к работе с любыми данными!

# 1. Здоровый образец (Нормальный GC-состав, нет мутаций)
normal_dna = "ATGCCGTAGTCGATAGCGATGCTAGTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

# 2. Образец с риском онкологии (Высокий GC + мутация в промоторе TERT)
# Здесь мы имитируем позицию 228 с заменой на T
cancer_risk_dna = "G" * 228 + "T" + "G" * 50 

# 3. Образец с Болезнью Гентингтона (Длинный повтор CAG)
huntington_dna = "ATGC" + "CAG" * 45 + "TAGC"

# ТЕСТИРОВАНИЕ:
print("Checking Normal DNA:")
universal_genomic_scanner(normal_dna)

print("\nChecking Cancer Risk DNA:")
universal_genomic_scanner(cancer_risk_dna)

print("\nChecking Huntington DNA:")
universal_genomic_scanner(huntington_dna)

Checking Normal DNA:
--- Starting Genomic Analysis Pipeline ---
[STEP 1] GC-Content: 51.39%
--- Analysis Complete ---

Checking Cancer Risk DNA:
--- Starting Genomic Analysis Pipeline ---
[STEP 1] GC-Content: 99.64%
[STEP 2] ALERT: TERT Promoter Mutation (C228T) detected! Immortality switch ON.
--- Analysis Complete ---

Checking Huntington DNA:
--- Starting Genomic Analysis Pipeline ---
[STEP 1] GC-Content: 65.73%
[STEP 3] PATHOGENIC REPEAT: Found 45 CAG repeats (Huntington Risk)
--- Analysis Complete ---
```
