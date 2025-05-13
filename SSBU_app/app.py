from shiny import App, render, ui
from shared import df
from app_ui import app_ui
from utils import check_hardy_weinberg, analyze_genotype_distribution
import pandas as pd

def server(input, output, session):
    @output
    @render.ui
    def page_ui():
        if input.page() == "Úvod":
            return ui.TagList(
                ui.h2("Vitajte v SSBU aplikacii"),
                ui.p("obsah")
            )
        elif input.page() == "Hardy-Weinberg":
            return ui.TagList(
                ui.h2("Hardy-Weinbergova rovnováha"),
                ui.output_table("observed_values_table"),
                ui.output_table("expected_values_table"),
                ui.output_table("chi2_test_table"),
                ui.output_table("hw_results_table")
            )
        elif input.page() == "Genotypy a predispozície":
            return ui.TagList(
                ui.h2("Genotypy a predispozície k hemochromatóze"),
                ui.p("Analýza zastúpenia genotypov a predispozícií k hereditárnej hemochromatóze."),
                ui.h3("Percentuálne zastúpenie genotypov"),
                ui.output_table("genotype_distribution_table"),
                ui.h3("Predispozícia k hereditárnej hemochromatóze"),
                ui.output_table("hemochromatosis_risk_table"),
                ui.h3("Počet pacientov podľa predispozície"),
                ui.output_text("patient_counts")
            )
        else:
            return ui.TagList(
                ui.h2("Očistený dataset"),
                ui.output_table("data_table")
            )

    @output
    @render.table
    def data_table():
        return df

    def generate_hw_table(df, value_extractor):
        results = []
        for column in ["HFE C187G (H63D) [HFE]", "HFE A193T (S65C) [HFE]", "HFE G845A (C282Y) [HFE]"]:
            hw_result = check_hardy_weinberg(df, column)
            results.append(value_extractor(hw_result, column) if hw_result else value_extractor(None, column))
        return pd.DataFrame(results)

    @output
    @render.table
    def hw_results_table():
        return generate_hw_table(df, lambda r, c: {
            "Mutácia": c,
            "Chi-štvorcový test": round(r["chi2"], 4) if r and r.get("chi2") is not None else "N/A",
            "p-hodnota": round(r["p_value"], 4) if r and r.get("p_value") is not None else "N/A",
            "Výsledok": "V rovnováhe" if r and r.get("p_value") is not None and r["p_value"] > 0.05 else "Nie je v rovnováhe" if r else "N/A"
        })
    
    @output
    @render.table
    def observed_values_table():
        return generate_hw_table(df, lambda r, c: {
            "Mutácia": c,
            "Normal": r["observed"]["normal"] if r else "N/A",
            "Heterozygot": r["observed"]["heterozygot"] if r else "N/A",
            "Mutant": r["observed"]["mutant"] if r else "N/A"
        })

    @output
    @render.table
    def expected_values_table():
        return generate_hw_table(df, lambda r, c: {
            "Mutácia": c,
            "Normal": round(r["expected"]["normal"], 2) if r else "N/A",
            "Heterozygot": round(r["expected"]["heterozygot"], 2) if r else "N/A",
            "Mutant": round(r["expected"]["mutant"], 2) if r else "N/A"
        })

    @output
    @render.table
    def chi2_test_table():
        return generate_hw_table(df, lambda r, c: {
            "Mutácia": c,
            "Chi-štvorcový test": round(r["chi2"], 4) if r else "N/A",
            "Stupne voľnosti": r["df"] if r else "N/A",
            "p-hodnota": round(r["p_value"], 4) if r and r["p_value"] else "N/A"
        })

    @output
    @render.table
    def genotype_distribution_table():
        return analyze_genotype_distribution(df)["genotype_percentages"]

    @output
    @render.table
    def hemochromatosis_risk_table():
        return analyze_genotype_distribution(df)["risk_percentages"]

    @output
    @render.text
    def patient_counts():
        results = analyze_genotype_distribution(df)
        total = len(df)
        carriers = results["carriers_count"]
        at_risk = results["at_risk_count"]
        return (
            f"Celkový počet pacientov: {total}\n"
            f"Počet prenášačov (heterozygoti): {carriers} ({carriers / total * 100:.2f}%)\n"
            f"Počet pacientov s genetickou predispozíciou: {at_risk} ({at_risk / total * 100:.2f}%)"
        )

app = App(app_ui, server)
