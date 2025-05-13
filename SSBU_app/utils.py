import pandas as pd
import numpy as np
from scipy.stats import chi2

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
        "Skupina": ["Prenášači", "Genetická predispozícia", "Celkový počet pacientov"],
        "Počet": [total_carriers, total_at_risk, total_patients],
        "Percento (%)": [round(total_carriers / total_patients * 100, 2),
                         round(total_at_risk / total_patients * 100, 2),
                         100.00]
    }
    predisposition_df = pd.DataFrame(predisposition_data)

    return {
        "genotype_percentages": genotype_df,
        "risk_percentages": risk_df,
        "predisposition_table": predisposition_df
    }

def prirad_kapitolu_mkch10(kod):
    if pd.isna(kod):
        return "Neznáma"
    kod = kod.strip().upper()
    if len(kod) >= 3 and kod[0].isalpha():
        prve_pismeno = kod[0]
        zvysok = kod[1:]
        try:
            cislo = int(zvysok[:2]) if len(zvysok) >= 2 and zvysok[:2].isdigit() else -1
        except ValueError:
            cislo = -1

        if prve_pismeno == 'A' and 0 <= cislo <= 99:
            return "Infekčné a parazitárne choroby (A00-B99)"
        elif prve_pismeno == 'C' and 0 <= cislo <= 97 or prve_pismeno == 'D' and 0 <= cislo <= 48:
            return "Nádory (C00-D48)"
        elif prve_pismeno == 'D' and 50 <= cislo <= 89:
            return "Choroby krvi a krvotvorných orgánov a niektoré poruchy imunity (D50-D89)"
        elif prve_pismeno == 'E' and 0 <= cislo <= 90:
            return "Endokrinné, nutričné a metabolické choroby (E00-E90)"
        elif prve_pismeno == 'F' and 0 <= cislo <= 99:
            return "Duševné poruchy a poruchy správania (F00-F99)"
        elif prve_pismeno == 'G' and 0 <= cislo <= 99:
            return "Choroby nervovej sústavy (G00-G99)"
        elif prve_pismeno == 'H' and 0 <= cislo <= 59:
            return "Choroby oka a adnexov (H00-H59)"
        elif prve_pismeno == 'H' and 60 <= cislo <= 95:
            return "Choroby ucha a mastoidného výbežku (H60-H95)"
        elif prve_pismeno == 'I' and 0 <= cislo <= 99:
            return "Choroby obehovej sústavy (I00-I99)"
        elif prve_pismeno == 'J' and 0 <= cislo <= 99:
            return "Choroby dýchacej sústavy (J00-J99)"
        elif prve_pismeno == 'K' and 0 <= cislo <= 93:
            return "Choroby tráviacej sústavy (K00-K93)"
        elif prve_pismeno == 'L' and 0 <= cislo <= 99:
            return "Choroby kože a podkožného tkaniva (L00-L99)"
        elif prve_pismeno == 'M' and 0 <= cislo <= 99:
            return "Choroby svalovej a kostrovej sústavy a spojivového tkaniva (M00-M99)"
        elif prve_pismeno == 'N' and 0 <= cislo <= 99:
            return "Choroby močovej a pohlavnej sústavy (N00-N99)"
        elif prve_pismeno == 'O' and 0 <= cislo <= 99:
            return "Tehotenstvo, pôrod a popôrodie (O00-O99)"
        elif prve_pismeno == 'P' and 0 <= cislo <= 96:
            return "Niektoré stavy vzniknuté v perinatálnom období (P00-P96)"
        elif prve_pismeno == 'Q' and 0 <= cislo <= 99:
            return "Vrodené chyby, deformity a chromozómové abnormality (Q00-Q99)"
        elif prve_pismeno == 'R' and 0 <= cislo <= 99:
            return "Príznaky, znaky a abnormálne klinické a laboratórne nálezy nezatriedené inde (R00-R99)"
        elif prve_pismeno == 'S' and 0 <= cislo <= 99 or prve_pismeno == 'T' and 0 <= cislo <= 98:
            return "Poranenia, otravy a niektoré iné následky vonkajších príčin (S00-T98)"
        elif prve_pismeno == 'V' and 0 <= cislo <= 99 or prve_pismeno == 'W' and 0 <= cislo <= 99 or prve_pismeno == 'X' and 0 <= cislo <= 99 or prve_pismeno == 'Y' and 0 <= cislo <= 99:
            return "Vonkajšie príčiny morbidity a mortality (V01-Y98)"
        elif prve_pismeno == 'Z' and 0 <= cislo <= 99:
            return "Faktory ovplyvňujúce zdravotný stav a kontakt so zdravotníckymi službami (Z00-Z99)"
        else:
            return "Iné nezatriedené"
    else:
        return "Neplatný kód"

def analyzuj_diagnozy(df, stlpec_mkch, stlpec_datum):
    df_analyza = df.copy()
    df_analyza['Rok_vyšetrenia'] = pd.to_datetime(df_analyza[stlpec_datum], format='%d.%m.%Y %H:%M', errors='coerce').dt.year
    df_analyza['Kapitola_MKCH10'] = df_analyza[stlpec_mkch].apply(prirad_kapitolu_mkch10)
    vyskyt_diagnoz = df_analyza.groupby(['Rok_vyšetrenia', 'Kapitola_MKCH10']).size().reset_index(name='Počet')
    return vyskyt_diagnoz