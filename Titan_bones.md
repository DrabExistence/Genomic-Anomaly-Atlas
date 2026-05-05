# Titan Bones Gene: LRP5 (D171V)
Description: The mutation causes the Wnt signaling pathway to be constantly active, resulting in extremely high bone density. People with this mutation are literally "unbreakable."

Detection Code:
Python
```python
def check_titan_bones(sequence, pos_in_seq):
    # Мутация D171V: замена аспартата на валин. 
    # На уровне ДНК это обычно замена A на T в определенном кодоне.
    base = sequence[pos_in_seq].upper()
    if base == 'T':
        return "🦾 TITAN BONES: High bone density mutation (LRP5) detected!"
    return "Normal bone density."
```
