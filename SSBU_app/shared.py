from pathlib import Path
import pandas as pd

app_dir = Path(__file__).parent
df = pd.read_csv(app_dir / "ocisteny_dataset.csv", sep=";", header=0, low_memory=False)