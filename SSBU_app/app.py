from shiny import App, render, ui
from shared import df
from app_ui import app_ui
from utils import check_hardy_weinberg
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
                ui.output_table("hw_table")
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

    @output
    @render.table
    def hw_table():
        results = {
            "Mutácia": [],
            "p-hodnota": [],
            "Frekvencia alely p": [],
            "Frekvencia alely q": [],
            "Pozorované normálne": [],
            "Pozorované heterozygot": [],
            "Pozorované mutant": [],
            "Očakávané normálne": [],
            "Očakávané heterozygot": [],
            "Očakávané mutant": []
        }

        for column in ["HFE C187G (H63D) [HFE]", "HFE A193T (S65C) [HFE]", "HFE G845A (C282Y) [HFE]"]:
            hw_result = check_hardy_weinberg(df, column)

            results["Mutácia"].append(column)

            if hw_result is None:
                results["p-hodnota"].append("N/A")
                results["Frekvencia alely p"].append("N/A")
                results["Frekvencia alely q"].append("N/A")
                for key in ["Pozorované normálne", "Pozorované heterozygot", "Pozorované mutant",
                            "Očakávané normálne", "Očakávané heterozygot", "Očakávané mutant"]:
                    results[key].append("N/A")
            else:
                results["p-hodnota"].append(
                    round(hw_result["p_value"], 4) if hw_result["p_value"] is not None else "N/A")
                results["Frekvencia alely p"].append(round(hw_result["allele_p"], 4))
                results["Frekvencia alely q"].append(round(hw_result["allele_q"], 4))

                results["Pozorované normálne"].append(hw_result["observed"]["normal"])
                results["Pozorované heterozygot"].append(hw_result["observed"]["heterozygot"])
                results["Pozorované mutant"].append(hw_result["observed"]["mutant"])

                results["Očakávané normálne"].append(round(hw_result["expected"]["normal"], 2))
                results["Očakávané heterozygot"].append(round(hw_result["expected"]["heterozygot"], 2))
                results["Očakávané mutant"].append(round(hw_result["expected"]["mutant"], 2))

        return pd.DataFrame(results)


app = App(app_ui, server)