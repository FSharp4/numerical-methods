import settings


def csv_read(file: str):
    with open(file, 'r') as csv:
        return csv.read()


def read_m19_bh():
    with open(f"{settings.PROJECT_ROOT}/m19_bh.csv", 'r') as file:
        return file.read()


if __name__ == "__main__":
    print("Test read:")
    print(read_m19_bh())
