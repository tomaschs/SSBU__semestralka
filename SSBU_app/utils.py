import pandas as pd
from scipy.stats import chi2_contingency
import numpy as np


def get_genotype_counts(df, column):
    counts = df[column].value_counts()
    n_normal = counts.get('normal', 0)              # dominantny homozygot p^2
    n_heterozygot = counts.get('heterozygot', 0)    # heterozygot 2pq
    n_mutant = counts.get('mutant', 0)              # recesivny homozygot q^2

    total = n_normal + n_heterozygot + n_mutant

    return n_normal, n_heterozygot, n_mutant, total


def perform_chi_square_test(observed, expected):
    # if min(expected) < 5:
    #     return {
    #         'chi2': None,
    #         'p_value': None,
    #         'df': None,
    #         'is_valid': False,
    #         'message': 'Expected counts too small for chi-square test (minimum count should be ≥ 5)'
    #     }

    chi2, p_value, dof, _ = chi2_contingency([observed, expected], correction=False)

    is_in_equilibrium = p_value > 0.05

    return {
        'chi2': chi2,
        'p_value': p_value,
        'df': dof,
        'is_valid': True,
        'message': 'V Hardy-Weinbergovej rovnováhe' if is_in_equilibrium else 'Nie je v Hardy-Weinbergovej rovnováhe'
    }


def check_hardy_weinberg(df, column):
    n_normal, n_heterozygot, n_mutant, total = get_genotype_counts(df, column)

    if n_normal + n_heterozygot + n_mutant < 2 or total == 0:
        return None

    # p + q = 1
    p = (2 * n_normal + n_heterozygot) / (2 * total)
    q = 1 - p

    expected = {
        'normal': (p ** 2) * total,  # p^2
        'heterozygot': 2 * p * q * total,  # 2pq
        'mutant': (q ** 2) * total  # q^2
    }

    observed = [n_normal, n_heterozygot, n_mutant]
    expected_values = [expected['normal'], expected['heterozygot'], expected['mutant']]

    chi_square_result = perform_chi_square_test(observed, expected_values)

    result = {
        'allele_p': p,
        'allele_q': q,
        'observed': dict(zip(['normal', 'heterozygot', 'mutant'], observed)),
        'expected': expected,
        'chi_square': chi_square_result
    }

    if chi_square_result['is_valid']:
        result['p_value'] = chi_square_result['p_value']
    else:
        result['p_value'] = None
        result['error'] = chi_square_result['message']

    return result


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
                      "C282Y/H63D zložený heterozygot", "C282Y/S65C zložený heterozygot",
                      "Prenášači", "Bez rizika"],
        "Počet": [n_mutant_c282y, n_mutant_h63d, n_mutant_s65c,
                  compound_heterozygotes, compound_c282y_s65c, total_carriers,
                  total_patients - total_carriers - total_at_risk],
        "Percento (%)": [
            n_mutant_c282y / total_patients * 100,
            n_mutant_h63d / total_patients * 100,
            n_mutant_s65c / total_patients * 100,
            compound_heterozygotes / total_patients * 100,
            compound_c282y_s65c / total_patients * 100,
            total_carriers / total_patients * 100,
            (total_patients - total_carriers - total_at_risk) / total_patients * 100
        ],
    }

    risk_df = pd.DataFrame(risk_data)
    risk_df["Percento (%)"] = risk_df["Percento (%)"].round(2)

    return {
        "genotype_percentages": genotype_df,
        "risk_percentages": risk_df,
        "carriers_count": total_carriers,
        "at_risk_count": total_at_risk
    }