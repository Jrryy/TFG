import csv
import random
from string import ascii_letters

with open("random_dataset.csv", "w") as f:
    csv_writer = csv.writer(f)

    headers = [" ", "class"] + list(
        {"gene_" + "".join(random.choices(ascii_letters, k=20)) for _ in range(20000)}
    )
    csv_writer.writerow(headers)
    length_genes = len(headers) - 2

    individuals = list(
        {"person_" + "".join(random.choices(ascii_letters, k=10)) for _ in range(100)}
    )
    length_individuals = len(individuals)

    for i in range(length_individuals):
        person = individuals[i]
        class_type = "AAA" if i < length_individuals/2 else "BBB"

        genes_expression = [random.uniform(-100, 100) for _ in range(length_genes)]

        csv_writer.writerow([person, class_type] + genes_expression)
