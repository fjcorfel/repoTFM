with open('data/mysplits.txt') as file:
    line = file.readline()

numbers = line.strip()[1:-1].split(', ')

with open('mysplits_column.txt', 'w') as output_file:
    for number in numbers:
        output_file.write(f"{number.strip()}\n")
        