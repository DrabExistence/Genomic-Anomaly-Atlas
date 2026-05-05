# Visualization - A Window into the Genome

Today we'll write a script that doesn't just output text, but also draws a graph. The graph will show how the GC composition changes along a DNA strand. To do this, we'll use the "sliding window" method.

How does "sliding window" work?
Imagine viewing a long sequence through a short frame (window) 50 letters long.

You calculate the GC % in this frame.

Move the frame one letter to the right.

Calculate again.

And so on until the end.

The result is a smooth curve that shows the peaks and valleys of genome density.

<img width="850" height="470" alt="image" src="https://github.com/user-attachments/assets/fbd6a0dc-8c43-4426-81ea-13d2bd15631d" />

# Python
```python
import matplotlib.pyplot as plt

def gc_content_analysis(sequence, window_size=50):
    gc_values = []
    
    # Скользящее окно
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        gc_count = window.count('G') + window.count('C')
        gc_percent = (gc_count / window_size) * 100
        gc_values.append(gc_percent)
        
    return gc_values

# Генерируем тестовую последовательность: 
# сначала обычная ДНК, потом CpG-островок (70% GC), потом снова обычная
test_dna = "AT" * 100 + "GC" * 50 + "GGC" * 30 + "AT" * 100

gc_data = gc_content_analysis(test_dna)

# Рисуем график
plt.figure(figsize=(10, 5))
plt.plot(gc_data, color='forestgreen', linewidth=2)
plt.axhline(y=41, color='red', linestyle='--', label='Genome Average (41%)')
plt.title("GC Content Profile: Identifying Promoter Zones")
plt.xlabel("Position in DNA (bp)")
plt.ylabel("GC %")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```
