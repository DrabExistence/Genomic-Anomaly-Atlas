# 🧬 Genomic Anomaly Atlas

**Genomic Anomaly Atlas** — это профессиональный конвейер (pipeline) для анализа генетических вариаций. Проект объединяет биоинформатическую обработку данных, поиск клинических соответствий в NCBI ClinVar и автоматическую генерацию PDF-отчетов.


---


## 🚀 Основные возможности v4.5

* **VCF Support:** Прямая обработка файлов Variant Call Format (VCF).
* **Exon Analysis:** Дифференциация изоформ FGFR2 (анализ сходства IIIb и IIIc).
* **Genomic Landscape:** Визуализация GC-состава («Генетический радар») для поиска аномалий.
* **ClinVar Integration:** Проверка патогенности вариантов через официальный API.
* **Auto-Reporting:** Создание PDF-отчетов с графиками.


---


## 🏛 Структура проекта

* `Genomic_Anomaly_Atlas.ipynb` — Основная интерактивная среда (Google Colab).
* `Genomic_Pipeline.py` — Универсальный скрипт для запуска в терминале или Docker.
* **modules/** — Исследовательские модули и база знаний (FGFR2, Huntington, ClinVar API).
* `test_data.vcf` — Тестовый образец геномных вариаций в формате VCF.
* `requirements.txt` — Список зависимостей для быстрой установки окружения.
* `sample_dna.txt` — (Опционально) Входной файл с сырой последовательностью ДНК.
* `rs_numbers.txt` — (Опционально) Список RS-идентификаторов для анализа.


---


## 🛠 Установка и запуск


### 1. Требования

Для работы необходим Python 3.8+ и следующие библиотеки:

```bash
pip install -r requirements.txt
```
---



### 2. Запуск через Терминал (CLI)


Скрипт поддерживает гибкую настройку через аргументы командной строки. Это основной способ интеграции в автоматизированные системы.


**Анализ VCF файла:**

```bash
python Genomic_Pipeline.py --vcf test_data.vcf --output Result.pdf
```
**Анализ по списку RS-номеров (из файла):**

```bash
python Genomic_Pipeline.py --input rs_numbers.txt
```
### 3. Работа в Google Colab

Вы можете запустить интерактивную версию проекта прямо в браузере. Она поддерживает визуализацию графиков прямо внутри блокнота:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DrabExistence/Genomic-Anomaly-Atlas/blob/main/Genomic_Anomaly_Atlas.ipynb)

**🐳 Docker (Изолированная среда)**
Использование Docker гарантирует наличие всех необходимых системных зависимостей и шрифтов:

```bash
docker build -t genomic-atlas .
docker run -v $(pwd):/app genomic-atlas
```
**👨‍🔬 Автор**
**AI-Bioinformatics Architect** — Разработка систем автоматизированной интерпретации геномных данных.

Проект создан для упрощения анализа редких генетических вариаций и наглядной визуализации геномного ландшафта.
