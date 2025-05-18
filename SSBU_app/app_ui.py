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
    .horizontal-buttons {
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
    .graf-kontajner {
        display: grid;
        grid-template_columns: repeat(2, 1fr);
        grid-gap: 20px;
    }
    .mutacia-sekcia {
        margin-bottom: 20px;
        vertical-align: top;
    }
    .graf-obrazok {
        margin: 10px 5px;
        display: block;
        max-width: 100%;
        height: auto;
    }

    /* Štýly pre menu (radio buttony) */
    .shiny-input-radiogroup .shiny-options-group {
        display: flex;
        flex-direction: column;
        width: 100%;
    }

    .shiny-input-radiogroup label {
        display: block;
        padding: 8px 12px;
        background-color: #f0f0f0;
        border: 1px solid #ccc;
        border-radius: 5px;
        cursor: pointer;
        transition: background-color 0.3s, color 0.3s;
        margin-bottom: 5px;
        text-align: center;
        width: 100%;
        box-sizing: border-box;
        margin-top: 0px;
    }

    .shiny-input-radiogroup label:hover,
    .shiny-input-radiogroup input[type="radio"]:checked + label {
        background-color: #007bff;
        color: white;
        border-color: #0056b3;
    }

    .shiny-input-radiogroup input[type="radio"] {
        display: none;
    }

    .shiny-input-radiogroup label a {
        text-decoration: none;
        color: inherit;
        display: block;
        width: 100%;
        box-sizing: border-box;
    }

""")

app_ui = ui.page_fluid(
    css,
    ui.layout_sidebar(
        ui.sidebar(
            ui.panel_title("HH-Predict"),
            ui.input_radio_buttons(
                "page", None,
                choices=["Úvod", "Data", "Hardy-Weinberg", "Genotypy a predispozície", "Analýza diagnóz", "MKCH-10", "Genotypy, demografia a diagnózy"],
                selected="Úvod",
                inline=False,
                width="100%"
            ),
            ui.output_ui("dynamic_content")
        ),
        ui.card(
            ui.output_ui("page_ui")
        )
    )
)