# LRP5: High Bone Mass (Titan Bones)

## 🧬 Biology
The **LRP5** gene is involved in the Wnt signaling pathway, which controls bone formation. 

## 🔍 Mutation: D171V (Aspartate to Valine)
A rare "gain-of-function" mutation at position 171. This makes the bone-building process nearly unstoppable, resulting in bones that are almost impossible to break.

## 💻 Logic
The pipeline scans the LRP5 sequence for a specific amino acid change (D->V) caused by a nucleotide substitution.

# Detection Code: Python
```python
def check_titan_bones(sequence, pos_in_seq):
    # Мутация D171V: замена аспартата на валин. 
    # На уровне ДНК это обычно замена A на T в определенном кодоне.
    base = sequence[pos_in_seq].upper()
    if base == 'T':
        return "🦾 TITAN BONES: High bone density mutation (LRP5) detected!"
    return "Normal bone density."
```
