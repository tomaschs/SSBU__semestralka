from shiny import ui

# CSS styles
css = ui.tags.style("""
    table.dataframe {
        width: 100%;
        border-collapse: collapse;
    }
    table.dataframe th, table.dataframe td {
        text-align: center;
        vertical-align: middle;
        border: 1px solid #ddd;
        padding: 8px;
    }
    table.dataframe th {
        background-color: #f2f2f2;
    }
""")

# Main UI definition
app_ui = ui.page_fluid(
    css,
    ui.layout_sidebar(
        ui.sidebar(
            ui.panel_title("Biomedical Data Visualization and Analysis"),
            ui.input_radio_buttons(
                "page", "Menu",
                choices=["Úvod", "Data"],
                selected="Úvod"
            ),
            ui.output_ui("dynamic_content")
        ),
        ui.card(
            ui.output_ui("page_ui")
        )
    )
)