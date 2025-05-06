
from shiny import ui
import data_init as data

app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.sidebar(
            ui.panel_title("CSV Data Viewer"),
            ui.input_text("search", "Search in data", placeholder="Type to search..."),
        ),
        ui.panel_main(
            ui.output_table("data_table")
        )
    )
)