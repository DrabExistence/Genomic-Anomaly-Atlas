Markdown
# Huntington's Disease: Trinucleotide Repeat Expansion

## 🧬 Theory
The **HTT** gene contains a sequence of three nucleotides (**CAG**) that repeat multiple times. 
* **Normal:** < 35 repeats.
* **Pathogenic:** > 36 repeats (Expansion).

## ⚠️ Consequence
The expanded polyglutamine (polyQ) tract causes the Huntingtin protein to misfold and aggregate, leading to progressive neuronal death in the brain.

## 🤖 AI Logic Task
Calculate the number of CAG repeats in a given sequence and determine the risk category.

# Python
```python
def check_huntington_risk(sequence):
# Считаем количество подряд идущих CAG
    import re
     repeats = re.findall(r'(?:CAG)+', sequence)
      if not repeats: return 0    
      max_repeats = len(max(repeats, key=len)) // 3    
     if max_repeats < 35:
        status = "Normal"
      elif 36 <= max_repeats <= 39:
        status = "Reduced Penetrance (Risk)"
      else:
        status = "Pathogenic"
        
    return max_repeats, status
    dna_sample = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" # Пример
    print(f"Repeats: {check_huntington_risk(dna_sample)}")
```
