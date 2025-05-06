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
        'is_in_equilibrium': is_in_equilibrium,
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