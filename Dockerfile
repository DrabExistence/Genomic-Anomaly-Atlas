Dockerfile
# Используем легкий образ Python
FROM python:3.9-slim

# Создаем папку для приложения
WORKDIR /app

# Копируем твои скрипты в контейнер
COPY . .

# Если нужны библиотеки (например, matplotlib)
RUN pip install matplotlib

# Команда для запуска пайплайна
CMD ["python", "Genomic_Pipeline.py"]
