import pandas as pd

data = pd.read_csv("data/ocisteny_dataset.csv", delimiter=";")
columns = data.columns.tolist()