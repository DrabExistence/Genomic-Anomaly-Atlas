# FGFR2 Isoforms Analysis: Epithelial vs Mesenchymal Switch

## 🧬 Theoretical Overview
This project explores the alternative splicing of the **FGFR2** (Fibroblast Growth Factor Receptor 2) gene. 
The switch between isoforms **IIIb** and **IIIc** is a fundamental mechanism in mammalian development.

### 🔍 Key Biological Insights
* **Isoform IIIb (Exon 8):** Expressed in epithelial cells. Essential for skin and organ development.
* **Isoform IIIc (Exon 9):** Expressed in mesenchymal cells. Crucial for bone and connective tissue formation.
* **Pathology:** Imbalances in this switch are linked to various cancers (Epithelial-Mesenchymal Transition) and genetic syndromes like Pfeiffer or Apert syndrome.

## 💻 Computational Analysis
I used an AI-assisted approach to compare the nucleotide sequences of Exon 8 and Exon 9.

### Similarity Script (Python)
```python
def calculate_similarity(seq1, seq2):
    # Приводим к одной длине для сравнения
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
    
    similarity = (matches / max(len(seq1), len(seq2))) * 100
    return round(similarity, 2)

# Пример последовательностей экзонов (упрощенно)
exon_8_IIIb = "ATTATAGTAGAGAGA" 
exon_9_IIIc = "ATTAGAGTAGAGACA"

result = calculate_similarity(exon_8_IIIb, exon_9_IIIc)
print(f"Сходство между изоформами IIIb и IIIc: {result}%")
