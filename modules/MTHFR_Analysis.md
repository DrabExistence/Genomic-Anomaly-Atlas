# Theory:
C677T polymorphism in the MTHFR gene.

Normal (C): The enzyme is 100% functional, and folic acid is absorbed.

Mutation (T): The enzyme is "weak" (functioning at 30-70%), which can lead to problems with utensils.

<img width="2000" height="2000" alt="image" src="https://github.com/user-attachments/assets/2e1393ea-9f3a-4ffe-b156-5f4978e12f85" />

# Python

Here the code checks not just the presence of a letter, but the genotype - after all, a person has two copies of each gene, from the father and from the mother

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
