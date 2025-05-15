from pathlib import Path
import pandas as pd

app_dir = Path(__file__).parent
df = pd.read_csv(app_dir / "ocisteny_dataset.csv", sep=";", header=0, low_memory=False)

# Cesta k súboru s MKCH-10 číselníkom
mkch10_file_path = app_dir / "Medzinarodna_klasifikacia_chorob_01062024.xls"
nazvy_harkov_mkch10 = ["A00-B99", "C00-D48", "D50-D90", "E00-E90", "F00-F99", "G00-G99", "H00-H59", "H60-H95", "I00-I99", "J00-J99", "K00-K93", "L00-L99", "M00-M99", "N00-N99", "O00-O99", "P00-P96", "Q00-Q99", "R00-R99", "S00-T98", "V01-Y98", "Z00-Z99", "U00-U99"]