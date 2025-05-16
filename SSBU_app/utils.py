import pandas as pd
import numpy as np
from scipy.stats import chi2
import logging

logging.basicConfig(level=logging.INFO)

def get_genotype_counts(df, column):
    counts = df[column].value_counts()
    n_normal = counts.get('normal', 0)
    n_heterozygot = counts.get('heterozygot', 0)
    n_mutant = counts.get('mutant', 0)
    total = n_normal + n_heterozygot + n_mutant
    return n_normal, n_heterozygot, n_mutant, total

def chi_square_test(observed, expected):
    observed = np.array(observed)
    expected = np.array(expected)

    valid_indices = expected > 0
    observed_valid = observed[valid_indices]
    expected_valid = expected[valid_indices]

    num_valid_categories = len(observed_valid)

    if num_valid_categories < 2:
        return 0, 1, 0 # Chi2 = 0, p=1, df=0 if less than 2 valid categories

    chi2_statistic = np.sum((observed_valid - expected_valid)**2 / expected_valid)
    degrees_of_freedom = 1 # Pre HWE s dvoma alelami sú vždy 1 stupeň voľnosti

    p_value = chi2.sf(chi2_statistic, degrees_of_freedom)

    print(f"Calculated Chi-squared statistic: {chi2_statistic}")
    print(f"Calculated Degrees of Freedom: {degrees_of_freedom}")
    print(f"Calculated p-value: {p_value}")

    return chi2_statistic, p_value, degrees_of_freedom


def check_hardy_weinberg(df, column):
    counts = df[column].value_counts()
    print(f"Value counts for column '{column}':\n{counts}")
    n_normal = counts.get('normal', 0)
    n_heterozygot = counts.get('heterozygot', 0)
    n_mutant = counts.get('mutant', 0)
    total = n_normal + n_heterozygot + n_mutant
    observed_counts = [n_normal, n_heterozygot, n_mutant]
    non_zero_observed = sum(1 for count in observed_counts if count > 0)

    print(f"n_normal: {n_normal}, n_heterozygot: {n_heterozygot}, n_mutant: {n_mutant}, total: {total}")
    print(f"Number of non-zero observed genotypes: {non_zero_observed}")

    reason_na = None
    allele_info = None
    expected = None
    chi2_results = None

    if total < 2:
        reason_na = "Príliš málo vzoriek"
        allele_info = {}
        expected = {}
        chi2_results = {}
    elif non_zero_observed < 2:
        reason_na = "Príliš málo pozorovaných genotypov"
        allele_info = {}
        expected = {}
        chi2_results = {}
    else:
        total_alleles = 2 * total
        normal_alleles = 2 * n_normal + n_heterozygot
        mutant_alleles = 2 * n_mutant + n_heterozygot

        p = normal_alleles / total_alleles
        q = mutant_alleles / total_alleles
        print(f"Allele frequencies for column '{column}': p = {p:.4f}, q = {q:.4f}")

        allele_info = {
            'total_alleles': total_alleles,
            'normal_alleles': normal_alleles,
            'mutant_alleles': mutant_alleles,
            'allele_p': p,
            'allele_q': q
        }

        expected = {
            'normal': (p ** 2) * total,
            'heterozygot': 2 * p * q * total,
            'mutant': (q ** 2) * total
        }
        expected_counts = [expected['normal'], expected['heterozygot'], expected['mutant']]
        non_zero_expected = sum(1 for count in expected_counts if count > 0)
        print(f"Expected genotype counts for column '{column}': {expected}")
        print(f"Number of non-zero expected genotypes: {non_zero_expected}")

        chi2_statistic, p_value, df_returned = chi_square_test(observed_counts, expected_counts)
        print(f"DOF returned from chi_square_test for column '{column}': {df_returned}")
        chi2_results = {
            'chi2': chi2_statistic,
            'p_value': p_value,
            'df': df_returned
        }

    return {
        'reason_na': reason_na,
        'allele_info': allele_info,
        'observed': dict(zip(['normal', 'heterozygot', 'mutant'], observed_counts)),
        'expected': expected,
        'chi2_results': chi2_results
    }

def analyze_genotype_distribution(df):
    c282y_col = "HFE G845A (C282Y) [HFE]"
    h63d_col = "HFE C187G (H63D) [HFE]"
    s65c_col = "HFE A193T (S65C) [HFE]"
    data = df.copy()
    total_patients = len(data)

    n_normal_c282y, n_heterozygot_c282y, n_mutant_c282y, _ = get_genotype_counts(data, c282y_col)
    n_normal_h63d, n_heterozygot_h63d, n_mutant_h63d, _ = get_genotype_counts(data, h63d_col)
    n_normal_s65c, n_heterozygot_s65c, n_mutant_s65c, _ = get_genotype_counts(data, s65c_col)

    genotype_percentages_data = {
        "Mutácia": ["C282Y", "H63D", "S65C"],
        "Normal (%)": [
            n_normal_c282y / total_patients * 100,
            n_normal_h63d / total_patients * 100,
            n_normal_s65c / total_patients * 100
        ],
        "Heterozygot (%)": [
            n_heterozygot_c282y / total_patients * 100,
            n_heterozygot_h63d / total_patients * 100,
            n_heterozygot_s65c / total_patients * 100
        ],
        "Mutant (%)": [
            n_mutant_c282y / total_patients * 100,
            n_mutant_h63d / total_patients * 100,
            n_mutant_s65c / total_patients * 100
        ]
    }
    genotype_df = pd.DataFrame(genotype_percentages_data)
    for col in ["Normal (%)", "Heterozygot (%)", "Mutant (%)"]:
        genotype_df[col] = genotype_df[col].round(2)

    compound_heterozygotes = data[(data[c282y_col] == "heterozygot") &
                                  (data[h63d_col] == "heterozygot")].shape[0]

    compound_c282y_s65c = data[(data[c282y_col] == "heterozygot") &
                                     (data[s65c_col] == "heterozygot")].shape[0]

    c282y_carriers = data[(data[c282y_col] == "heterozygot") &
                          (data[h63d_col] != "heterozygot") &
                                (data[s65c_col] != "heterozygot")].shape[0]

    h63d_carriers = data[(data[h63d_col] == "heterozygot") &
                               (data[c282y_col] != "heterozygot") &
                                (data[s65c_col] != "heterozygot")].shape[0]

    s65c_carriers = data[(data[s65c_col] == "heterozygot") &
                               (data[c282y_col] != "heterozygot") &
                               (data[h63d_col] != "heterozygot")].shape[0]

    total_carriers = c282y_carriers + h63d_carriers + s65c_carriers
    total_at_risk = n_mutant_c282y + n_mutant_h63d + compound_heterozygotes + compound_c282y_s65c
    total_non_affected = total_patients - total_carriers - total_at_risk  # Calculate non-affected
    
    risk_data = {
        "Kategória": ["C282Y homozygot", "H63D homozygot", "S65C homozygot",
                      "C282Y/H63D zložený heterozygot", "C282Y/S65C zložený heterozygot"],
        "Percento (%)": [
            n_mutant_c282y / total_patients * 100,
            n_mutant_h63d / total_patients * 100,
            n_mutant_s65c / total_patients * 100,
            compound_heterozygotes / total_patients * 100,
            compound_c282y_s65c / total_patients * 100
        ],
    }
    risk_df = pd.DataFrame(risk_data)
    risk_df["Percento (%)"] = risk_df["Percento (%)"].round(2)

    predisposition_data = {
        "Skupina": ["Prenášači", "Genetická predispozícia", "Nepostihnutí", "Celkový počet pacientov"], # Add "Nepostihnutí"
        "Počet": [total_carriers, total_at_risk, total_non_affected, total_patients], # Add total_non_affected
        "Percento (%)": [round(total_carriers / total_patients * 100, 2),
                         round(total_at_risk / total_patients * 100, 2),
                         round(total_non_affected / total_patients * 100, 2), # Calculate non-affected %
                         100.00]
    }
    predisposition_df = pd.DataFrame(predisposition_data)

    return {
        "genotype_percentages": genotype_df,
        "risk_percentages": risk_df,
        "predisposition_table": predisposition_df
    }

def prirad_kapitolu_mkch10(kod, mkch10_data):
    """
    Priradí ku kódu MKCH-10 názov kategórie/podkategórie na základe dát z Excelu.

    Args:
        kod (str): Kód diagnózy MKCH-10.
        mkch10_data (dict):  Slovník, kde kľúče sú názvy hárkov Excelu
                             a hodnoty sú DataFrame s dátami MKCH-10.

    Returns:
        str: Názov kategórie/podkategórie alebo "Neznáma" ak sa kód nenájde.
    """
    if pd.isna(kod):
        return "Neznáma"
    kod = kod.strip().upper()

    for sheet_name, df in mkch10_data.items():
        kod_col = 'Kód diagnózy' if 'Kód diagnózy' in df.columns else 'Kód' if 'Kód' in df.columns else None
        nazov_col = 'Nazov' if 'Nazov' in df.columns else 'Názov' if 'Názov' in df.columns else None

        if kod_col and nazov_col:
            if kod in df[kod_col].values:
                return df.loc[df[kod_col] == kod, nazov_col].iloc[0]  # Vráti názov
    return "Neznáma"  # Ak sa kód nenájde v žiadnom hárku

def analyzuj_diagnozy(df, stlpec_mkch, stlpec_datum, mkch10_data):
    """
    Analyzuje výskyt diagnóz podľa MKCH-10 kódu a roka vyšetrenia,
    vracia aj percentuálne zastúpenie a informácie pre zvýraznenie chybných kódov.

    Args:
        df (pd.DataFrame): DataFrame s dátami pacientov.
        stlpec_mkch (str): Názov stĺpca s MKCH-10 kódom.
        stlpec_datum (str): Názov stĺpca s dátumom vyšetrenia.
        mkch10_data (dict):  Slovník, kde kľúče sú názvy hárkov Excelu
                             a hodnoty sú DataFrame s dátami MKCH-10.

    Returns:
        pd.DataFrame: DataFrame s počtom a percentuálnym zastúpením
                      jednotlivých MKCH-10 kódov pre každý rok vyšetrenia,
                      a zoznam chybných kódov.
    """

    df_analyza = df.copy()
    df_analyza['Rok_vyšetrenia'] = pd.to_datetime(df_analyza[stlpec_datum], format='%d.%m.%Y %H:%M', errors='coerce').dt.year
    df_analyza['Nazov_MKCH10'] = df_analyza[stlpec_mkch].apply(lambda x: prirad_kapitolu_mkch10(x, mkch10_data))

    # Počítame výskyty
    vyskyt_diagnoz = df_analyza.groupby(['Rok_vyšetrenia', stlpec_mkch, 'Nazov_MKCH10']).size().reset_index(name='Pocet')

    # Počítame percentá
    total_vyskytov = len(df_analyza)
    vyskyt_diagnoz['Percento'] = ((vyskyt_diagnoz['Pocet'] / total_vyskytov) * 100).round(2)

    # Hľadanie chybných kódov
    platne_kody = set()
    for sheet_name, sheet_df in mkch10_data.items():
        kod_col = 'Kód diagnózy' if 'Kód diagnózy' in sheet_df.columns else 'Kód' if 'Kód' in sheet_df.columns else None
        if kod_col:
            platne_kody.update(sheet_df[kod_col].astype(str).str.upper().tolist())

    vyskyt_diagnoz['Je_chybny'] = vyskyt_diagnoz[stlpec_mkch].astype(str).str.upper().apply(lambda x: x not in platne_kody)
    chybne_kody = vyskyt_diagnoz[vyskyt_diagnoz['Je_chybny'] == True][[stlpec_mkch]].drop_duplicates()  # Získame len unikátne chybné kódy

    return {
        'vyskyt_diagnoz': vyskyt_diagnoz,
        'celkovy_pocet': total_vyskytov,
        'chybne_kody': chybne_kody
    }

def nacitaj_mkch10_ciselnik(cesta_k_suboru, nazvy_harkov):
    try:
        all_sheets_dict = pd.read_excel(cesta_k_suboru, sheet_name=None) # Načíta všetky hárky
        logging.info(f"Úspešne načítaných {len(all_sheets_dict)} hárkov z číselníka MKCH-10.")
        return all_sheets_dict
    except FileNotFoundError:
        logging.error(f"Súbor '{cesta_k_suboru}' nebol nájdený.")
        return None
    except Exception as e:
        logging.error(f"Chyba pri načítaní číselníka MKCH-10: {e}")
        return None
    

