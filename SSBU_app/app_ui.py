from shiny import ui

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
        background-color: aliceblue;
    }
    #mkch10_hladaj {
        height: 30px;
        padding: 5px;
        font-size: 0.9em;
        margin-bottom: 5px;
    }
    .horizontal-buttons { /* Nový štýl pre tlačidlá vedľa seba */
        display: flex;
        flex-direction: row;
        align-items: center;
    }
    .horizontal-buttons > * {
        margin-right: 5px;
    }
    .horizontal-buttons > *:last-child {
        margin-right: 0;
    }
    #mkch10_hladaj_button, #mkch10_reset_search {
        padding: 5px 10px;
        font-size: 0.8em;
        height: 30px;
    }
    .horizontal-layout {
        display: flex;
        flex-direction: row;
        align-items: start;
        justify-content: space-between;
        margin-bottom: 10px;
    }
    .vertical-layout {
        display: flex;
        flex-direction: column;
        align-items: start;
    }
    .vertical-layout > * {
        margin-bottom: 5px;
    }
    .action-button-sm {
        padding: 5px 10px;
        font-size: 0.8em;
    }
    .action-button-nav {
        padding: 5px 10px;
        font-size: 0.9em;
        margin-bottom: 5px;
        cursor: pointer;
    }
""")

app_ui = ui.page_fluid(
    css,
    ui.layout_sidebar(
        ui.sidebar(
            ui.panel_title("SSBU"),
            ui.input_radio_buttons(
                "page", "Menu",
                choices=["Úvod", "Data", "Hardy-Weinberg", "Genotypy a predispozície", "Analýza diagnóz", "MKCH-10"], # Pridali sme "Analýza diagnóz"
                selected="Úvod"
            ),
            ui.output_ui("dynamic_content")
        ),
        ui.card(
            ui.output_ui("page_ui")
        )
    )
)