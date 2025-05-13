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

    if total < 2 or non_zero_observed < 2:
        print(f"Insufficient data for Hardy-Weinberg test for column '{column}'. Returning None.")
        return None

    p = (2 * n_normal + n_heterozygot) / (2 * total)
    q = 1 - p
    print(f"Allele frequencies for column '{column}': p = {p:.4f}, q = {q:.4f}")

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

    return {
        'allele_p': p,
        'allele_q': q,
        'observed': dict(zip(['normal', 'heterozygot', 'mutant'], observed_counts)),
        'expected': expected,
        'chi2': chi2_statistic,
        'p_value': p_value,
        'df': df_returned
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