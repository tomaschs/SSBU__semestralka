from shiny import App, render, ui
from shared import df
from app_ui import app_ui

def server(input, output, session):
    @output
    @render.ui
    def page_ui():
        if input.page() == "Úvod":
            return ui.TagList(
                ui.h2("Vitajte v SSBU aplikacii"),
                ui.p("obsah")
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

app = App(app_ui, server)