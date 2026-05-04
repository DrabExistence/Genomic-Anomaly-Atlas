<img width="2000" height="2000" alt="image" src="https://github.com/user-attachments/assets/2e1393ea-9f3a-4ffe-b156-5f4978e12f85" />

# Python
```python
 def analyze_mthfr(allele1, allele2):
  genotype = (allele1 + allele2).upper()    
    results = {
        "CC": "Normal (Wild type). Enzyme activity is 100%.",
        "CT": "Heterozygous mutation. Enzyme activity reduced by ~35%.",
        "TT": "Homozygous mutation. Enzyme activity reduced by ~70%. Risk of high homocysteine.",
        "TC": "Heterozygous mutation. Enzyme activity reduced by ~35%."
    }
    return results.get(genotype, "Unknown genotype. Please use C or T.") # Пример: проверка генотипа клиента
  print(f"Result for genotype CT: {analyze_mthfr('C', 'T')}")
```
