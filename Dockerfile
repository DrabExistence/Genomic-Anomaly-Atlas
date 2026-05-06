# Используем легкую версию Python
FROM python:3.9-slim

# Устанавливаем системные зависимости (шрифты и библиотеки для графики)
RUN apt-get update && apt-get install -y \
    fonts-dejavu-core \
    libfreetype6-dev \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Создаем рабочую директорию внутри контейнера
WORKDIR /app

# Копируем файл зависимостей
COPY requirements.txt .

# Устанавливаем библиотеки Python
RUN pip install --no-cache-dir -r requirements.txt

# Копируем все файлы проекта в контейнер
COPY . .

# Команда по умолчанию при запуске контейнера
# Она запустит твой основной пайплайн
CMD ["python", "Genomic_Pipeline.py"]
