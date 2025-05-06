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
            "p-hodnota": []
        }
        for column in ["HFE C187G (H63D) [HFE]", "HFE A193T (S65C) [HFE]", "HFE G845A (C282Y) [HFE]"]:
            p_value = check_hardy_weinberg(df, column)
            results["Mutácia"].append(column)
            results["p-hodnota"].append(p_value if p_value is not None else "N/A")
        return pd.DataFrame(results)

app = App(app_ui, server)