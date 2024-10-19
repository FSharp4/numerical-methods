import csv

FILE = 5

if __name__ == "__main__":
    with open(f"../circuits/circuit{FILE}.csv", newline='', encoding='utf-8-sig') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in csv_reader:
            print(row)
