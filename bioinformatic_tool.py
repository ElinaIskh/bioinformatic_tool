def run_dna_rna_tools(*args: str):
    """
    Takes 1 or more sequences, checks if is it nucleic acid and
    does one of the tools: is_nucleic_acid, reverse, transcribe,
    complement, reverse_complement

    Arguments:
    sequences: list[str]
    tool: str

    Returns bool/str
    Raises exception if:
        Number of argumets is less then 2
        Sequence is not an nucleic acid
        Tool is unknown
    """
    # Проверяем, что аргументов минимум 2
    if len(args) < 2:
        raise ValueError(
            "Введите хотя бы одну последовательность и одну процедуру")

    sequences = args[:-1]
    tool = args[-1]

    # словарь тулов
    tools = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    # Проверка, что процедура существует
    if tool not in tools:
        raise ValueError("Неизвестная процедура")

    # Проверка, что последовательность - нуклеиновая кислота
    for seq in sequences:
        if tool != "is_nucleic_acid" and not is_nucleic_acid(seq):
            raise ValueError("Введена некорректная последовательность")

    # Выполнение функции
    results = [tools[tool](seq) for seq in sequences]

    # Условие, чтобы при подаче 1 последовательности возвращалась строка
    if len(results) == 1:
        return results[0]
    return results
