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

    if min(expected_values) < 5:
        return {
            'p_value': None,
            'allele_p': p,
            'allele_q': q,
            'observed': dict(zip(['normal', 'heterozygot', 'mutant'], observed)),
            'expected': expected,
            'error': 'Expected counts too small for chi-square test'
        }

    # chi-kvadrat test
    chi2, p_value, _, _ = chi2_contingency([observed, expected_values], correction=False)

    return {
        'p_value': p_value,
        'allele_p': p,
        'allele_q': q,
        'observed': dict(zip(['normal', 'heterozygot', 'mutant'], observed)),
        'expected': expected
    }