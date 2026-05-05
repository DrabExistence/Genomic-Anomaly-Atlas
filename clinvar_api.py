!pip install tqdm
!pip install biopython

from Bio import Entrez
import xml.etree.ElementTree as ET
from urllib.error import URLError, HTTPError
from tqdm.auto import tqdm

def get_clinvar_significance(rs_id):
    """
    Принимает rs-номер (например, 'rs1137282'),
    обращается к NCBI ClinVar и возвращает Clinical Significance.
    """
    Entrez.email = "example@gmail.com"  # Укажи любую почту

    try:
        # Шаг 1: Ищем ID в базе ClinVar по rs-номеру
        search_query = f"{rs_id}[Variant ID]"
        handle = Entrez.esearch(db="clinvar", term=search_query)
        search_results = Entrez.read(handle)
        handle.close()
    except (URLError, HTTPError) as e:
        return f"Network error during search: {e}"

    if not search_results['IdList']:
        return "No ClinVar record found."

    try:
        # Шаг 2: Получаем детали записи
        clinvar_id = search_results['IdList'][0]
        handle = Entrez.esummary(db="clinvar", id=clinvar_id)
        summary = Entrez.read(handle, validate=False)
        handle.close()
    except (URLError, HTTPError) as e:
        return f"Network error during summary retrieval: {e}"

    # Отладка: выводим структуру summary
    # print("Debug: Summary structure for rs_id", rs_id, ":\n", summary) # Закомментируем отладочный вывод для чистоты

    # Извлекаем клиническую значимость
    try:
        # В структуре ClinVar значение лежит в поле 'germline_classification' -> 'description'
        significance = summary['DocumentSummarySet']['DocumentSummary'][0]['germline_classification']['description']
        return significance
    except (KeyError, IndexError):
        return "Significance data unavailable."

# --- Оптимизированная внутренняя функция для пакетной обработки ---
def _get_clinvar_significance_batch_internal(rs_ids_list):
    """
    Внутренняя оптимизированная функция для пакетного получения клинической значимости.
    """
    Entrez.email = "example@gmail.com" # Убедиться, что email установлен
    final_results = {}
    clinvar_ids_to_fetch = []
    clinvar_id_map = {} # Сопоставление ClinVar ID с исходным rs_id

    # Шаг 1: Выполняем esearch для каждого rs_id, чтобы получить ClinVar ID
    for rs_id in tqdm(rs_ids_list, desc="Поиск ClinVar ID"):
        try:
            search_query = f"{rs_id}[Variant ID]"
            handle = Entrez.esearch(db="clinvar", term=search_query)
            search_results = Entrez.read(handle)
            handle.close()
        except (URLError, HTTPError) as e:
            final_results[rs_id] = f"Сетевая ошибка при поиске {rs_id}: {e}"
            continue

        if search_results['IdList']:
            clinvar_id = search_results['IdList'][0]
            clinvar_ids_to_fetch.append(clinvar_id)
            clinvar_id_map[clinvar_id] = rs_id
        else:
            final_results[rs_id] = "Запись ClinVar не найдена."

    # Шаг 2: Выполняем один запрос esummary для всех собранных ClinVar ID
    if clinvar_ids_to_fetch:
        id_list_str = ",".join(clinvar_ids_to_fetch)
        try:
            handle = Entrez.esummary(db="clinvar", id=id_list_str, validate=False)
            summary_list = Entrez.read(handle, validate=False)
            handle.close()
        except (URLError, HTTPError) as e:
            for clinvar_id in clinvar_ids_to_fetch:
                rs_id = clinvar_id_map[clinvar_id]
                final_results[rs_id] = f"Сетевая ошибка при получении сводки для {rs_id}: {e}"
            return final_results

        for doc_summary in summary_list['DocumentSummarySet']['DocumentSummary']:
            # Используем 'variation_set_id' как ключ для получения ClinVar ID из пакетного esummary
            clinvar_uid = doc_summary['variation_set_id']
            original_rs_id = clinvar_id_map.get(clinvar_uid)
            if original_rs_id:
                try:
                    significance = doc_summary['germline_classification']['description']
                    final_results[original_rs_id] = significance
                except (KeyError, IndexError):
                    final_results[original_rs_id] = "Данные о значимости недоступны."

    return final_results

# --- Оболочка для сохранения существующего имени функции ---
def get_clinvar_significance_for_list(rs_ids_list):
    """
    Принимает список rs-номеров и возвращает словарь с Clinical Significance для каждого,
    используя оптимизированные пакетные запросы.
    """
    return _get_clinvar_significance_batch_internal(rs_ids_list)


# Пример: проверка с одним rs-номером
print(f"ClinVar Status for rs1137282: {get_clinvar_significance('rs1137282')}")

# Пример: проверка со списком rs-номеров
rs_numbers_to_check = ['rs1137282', 'rs7069102', 'rs1234567'] # rs7069102 ранее не находился, rs1234567 - вымышленный для теста
list_results = get_clinvar_significance_for_list(rs_numbers_to_check)
print("\nClinVar Status for list of rs-numbers:")
for rs_id, significance in list_results.items():
    print(f"  {rs_id}: {significance}")
