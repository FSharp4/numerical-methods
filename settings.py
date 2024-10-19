import os
from pathlib import Path

__project_path__ = Path(os.path.realpath(__file__)).parent

PROJECT_ROOT = __project_path__.__str__()

if __name__ == "__main__":
    print(PROJECT_ROOT)