# Theory
Previously, 98% of our DNA was called "junk," but we now know that it contains switches (promoters and enhancers).

New Anomaly: Mutations in the TERT Gene Promoter (Cancer's "Immortality Button")

Theory: The TERT gene is responsible for completing the ends of chromosomes (telomeres). In normal cells, it is turned off, which is why we age. But in cancer cells, a "breakdown" occurs in the promoter (the region of DNA before the gene).

The anomaly: A simple change of one letter in the "junk" region causes the gene to turn on at full power. The cell stops aging and begins to divide endlessly, turning into a tumor.

For the Architect: This is an example of how a mutation outside the gene can be more dangerous than a mutation within the gene.

These mutations don't just change a letter. They create a new landing pad for proteins (ETS transcription factors).
Imagine a factory (gene) with a fence. A mutation is a hole in the fence through which "illegal" workers (proteins) begin to enter and start production (telomerase) when the factory should be closed.

# Python
```python
def analyze_tert_promoter(sequence, offset=0):
    """
    Сканирует промотор TERT на мутации C228T и C250T.
    sequence: строка ДНК (промоторная зона)
    offset: смещение, если последовательность начинается не с начала референса
    """
    # Горячие точки (позиции в референсной последовательности промотора)
    # C228T и C250T — это замены Цитозина на Тимин
    hotspots = {
        228: {'ref': 'C', 'alt': 'T', 'name': 'C228T'},
        250: {'ref': 'C', 'alt': 'T', 'name': 'C250T'}
    }
    
    found_mutations = []
    
    for pos, info in hotspots.items():
        rel_pos = pos - offset
        if rel_pos < len(sequence):
            actual_base = sequence[rel_pos].upper()
            if actual_base == info['alt']:
                found_mutations.append(info['name'])
            elif actual_base != info['ref']:
                print(f"Warning: Unexpected base at position {pos}: {actual_base}")

    # Интерпретация результатов
    if found_mutations:
        status = "⚠️ ACTIVATED (Oncogenic Risk)"
        desc = f"The 'immortality switch' is ON. Found: {', '.join(found_mutations)}. This mutation creates a new binding site for ETS transcription factors, overexpressing Telomerase."
    else:
        status = "✅ NORMAL"
        desc = "No hot-spot TERT promoter mutations detected. Telomerase expression remains regulated."

    return {
        "Status": status,
        "Details": desc,
        "Detected": found_mutations
    }

# --- Пример работы ---
# Представим кусок ДНК пациента с мутацией в 228 позиции
sample_dna = "G" * 228 + "T" + "G" * 21  # Упрощенно: в 228 позиции стоит T
report = analyze_tert_promoter(sample_dna)

print(f"TERT Analysis: {report['Status']}")
print(f"Clinical Note: {report['Details']}")
```
