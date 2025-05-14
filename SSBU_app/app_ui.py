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
    #mkch10_hladaj_button, #mkch10_zrusit_filter {
        padding: 5px 10px; /* Zmenšenie vnútorného priestoru */
        font-size: 0.8em; /* Zmenšenie písma */
    }
    .horizontal-layout {
        display: flex;
        flex-direction: row;
        align-items: center; /* Voliteľné: zarovnanie prvkov vertikálne na stred */
    }
    .horizontal-layout > * {
        margin-right: 5px; /* Voliteľné: pridanie medzier medzi prvkami */
    }
    .horizontal-layout > *:last-child {
        margin-right: 0; /* Odstránenie medzery za posledným prvkom */
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