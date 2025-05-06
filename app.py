
from shiny import App, render, reactive
import data_init as data
import utils

def server(input, output, session):
    @reactive.Calc
    def filtered_df():
        return utils.filter_data(data.data, input.search())

    @output
    @render.table
    def data_table():
        return filtered_df()

import app_ui
app = App(app_ui.app_ui, server)