import logging
from shiny import App, render, ui, reactive
from shiny.ui import tags
import pandas as pd
from shared import df, mkch10_file_path, nazvy_harkov_mkch10
from app_ui import app_ui
from utils import check_hardy_weinberg, analyze_genotype_distribution, analyzuj_diagnozy, prirad_kapitolu_mkch10, \
    nacitaj_mkch10_ciselnik, chi_square_test
from plotnine import ggplot, aes, geom_bar, facet_wrap, geom_boxplot, geom_point, theme_minimal, labs, position_dodge, scale_fill_manual
from scipy.stats import chi2_contingency
import io
import base64
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.pyplot as plt
from utils import check_hardy_weinberg, analyze_genotype_distribution, analyzuj_diagnozy, prirad_kapitolu_mkch10, nacitaj_mkch10_ciselnik, format_p_value

mkch10_data = nacitaj_mkch10_ciselnik(str(mkch10_file_path), nazvy_harkov_mkch10)


def server(input, output, session):
    sheet_names = reactive.Value(list(mkch10_data.keys()) if mkch10_data else [])
    current_sheet_index = reactive.Value(0)
    filtered_mkch10_data = reactive.Value(None)
    search_performed = reactive.Value(False)
    filtered_data = reactive.Value(df)
    diagnozy = sorted(df["diagnoza MKCH-10"].dropna().unique())
    diagnozy.insert(0, "Všetky")

    def priprav_data_pre_grafy(df, vek_od=None, vek_do=None, pohlavie=None, diagnoza=None):
        data = df.copy()
        print(f"[DEBUG] Pôvodný počet riadkov: {len(data)}")

        if vek_od is not None:
            data = data[data['vek'] >= vek_od]
            print(f"[DEBUG] Po filtrovaní vek_od ({vek_od}): {len(data)} riadkov")

        if vek_do is not None:
            data = data[data['vek'] <= vek_do]
            print(f"[DEBUG] Po filtrovaní vek_do ({vek_do}): {len(data)} riadkov")

        mapovanie = {
            "Žena": "F",
            "Muž": "M",
            "Všetky": None
        }
        pohlavie_filter = mapovanie.get(pohlavie, None)

        if pohlavie_filter:
            data = data[data['pohlavie'] == pohlavie_filter]
            print(f"[DEBUG] Po filtrovaní pohlavie ({pohlavie}): {len(data)} riadkov")

        if diagnoza and diagnoza.strip() != "":
            data["diagnoza_ano_nie"] = data["diagnoza MKCH-10"].apply(
                lambda x: "Áno" if diagnoza in str(x) else "Nie"
            )
            print(f"[DEBUG] Pridaný stĺpec diagnoza_ano_nie pre diagnózu '{diagnoza}'")
        else:
            data["diagnoza_ano_nie"] = "Všetci"
            print("[DEBUG] Diagnóza nie je špecifikovaná, nastavené na 'Všetci'")

        print(f"[DEBUG] Celkový počet riadkov po filtrovaní: {len(data)}")
        return data

    def generuj_grafy(df, pohlavie=None, diagnoza=None, vek_od=None, vek_do=None):
        grafy_html = []
        mutacie = ["HFE G845A (C282Y) [HFE]", "HFE C187G (H63D) [HFE]", "HFE A193T (S65C) [HFE]"]

        df_copy = df.copy()

        if diagnoza and diagnoza.strip() != "" and diagnoza != "Všetky":
            df_copy["diagnoza_ano_nie"] = df_copy["diagnoza MKCH-10"].apply(
                lambda x: "Áno" if diagnoza in str(x) else "Nie"
            )
            fill_col = "diagnoza_ano_nie"
            fill_legend = f"Diagnóza {diagnoza}"
        else:
            df_copy["diagnoza_ano_nie"] = "Všetci"
            fill_col = None
            fill_legend = None

        title_filters = []
        if pohlavie and pohlavie != "Všetky":
            title_filters.append(f"Pohlavie: {pohlavie}")
        if diagnoza and diagnoza.strip() != "" and diagnoza != "Všetky":
            title_filters.append(f"Diagnóza: {diagnoza}")
        if vek_od is not None or vek_do is not None:
            vek_text = "Vek"
            if vek_od is not None:
                vek_text += f" od {vek_od}"
            if vek_do is not None:
                vek_text += f" do {vek_do}"
            title_filters.append(vek_text)

        filter_popis = ", ".join(title_filters) if title_filters else "Všetky dáta"

        def graf_na_html(plot_obj):
            fig = plot_obj.draw()
            buf = io.BytesIO()
            canvas = FigureCanvas(fig)
            canvas.print_png(buf)
            data = base64.b64encode(buf.getvalue()).decode("utf-8")
            plt.close(fig)
            return f'<img src="data:image/png;base64,{data}" class:"graf-obrazok" />'

        for mutacia in mutacie:
            sekcia_html = f'<div class="mutacia-sekcia" >'
            sekcia_html += f"<h2>{mutacia}</h2>"
            sekcia_html += f"<p>Filtre: {filter_popis}</p>"

            # 1. distribúcia genotypov – bez fill podla diagnozy; len čistý count
            p1 = ggplot(df_copy, aes(x=mutacia)) + geom_bar() \
                 + labs(title=f"{mutacia} – Distribúcia genotypov", x="Genotyp", y="Počet pacientov") \
                 + theme_minimal()
            sekcia_html += graf_na_html(p1)

            # 2. vek
            p2 = (ggplot(df_copy, aes(x=mutacia, y="vek"))
                  + geom_boxplot()
                  + labs(title=f"{mutacia} – Vek podľa genotypu", x="Genotyp", y="Vek")
                  + theme_minimal())
            sekcia_html += graf_na_html(p2)

            # 3. podla pohlavia
            p3 = (ggplot(df_copy, aes(x=mutacia, fill="pohlavie"))
                  + geom_bar(position=position_dodge())
                  + labs(title=f"{mutacia} – Genotypy podľa pohlavia", x="Genotyp", y="Počet pacientov",
                         fill="Pohlavie")
                  + theme_minimal())
            sekcia_html += graf_na_html(p3)

            # 4. podla diagnozy – len ak je diagnoza zvolena
            if fill_col:
                p4 = (ggplot(df_copy, aes(x=mutacia, fill=fill_col))
                      + geom_bar(position=position_dodge())
                      + labs(title=f"{mutacia} – Genotypy podľa diagnózy", x="Genotyp", y="Počet pacientov",
                             fill=fill_legend)
                      + theme_minimal())
                sekcia_html += graf_na_html(p4)
            sekcia_html += '</div>'
            grafy_html.append(sekcia_html)

        final_html = '<div class="graf-kontajner" >' + ''.join(grafy_html) + '</div>'
        return ui.HTML(final_html)

    def vykonaj_chi_kvadrat_testy(df, diagnoza=None, vek_od=None, vek_do=None):
        """
        Vykoná Chí-kvadrát testy na zistenie asociácie medzi HFE mutáciami
        a pohlavím/diagnózou/vekom, používajúc funkciu chi_square_test z utils.py.


        Args:
            df (pd.DataFrame): DataFrame s dátami pacientov.
            diagnoza (str, voliteľné): Konkrétna diagnóza na testovanie.
                                       Ak je None alebo prázdny, test sa vynechá.
            vek_od (int, voliteľné): Minimálny vek pre filtrovanie.
            vek_do (int, voliteľné): Maximálny vek pre filtrovanie.


        Returns:
            ui.TagList: HTML obsah s výsledkami testov.
        """

        vysledky = []
        mutacie = ["HFE G845A (C282Y) [HFE]", "HFE C187G (H63D) [HFE]", "HFE A193T (S65C) [HFE]"]
        pecen_kody = ['K76.0', 'K75.9']

        pecen_maska = df["diagnoza MKCH-10"].apply(
            lambda x: any(kod in str(x) for kod in pecen_kody)
        )

        for mutacia in mutacie:
            popis = f"{mutacia}"
            podmienky = []

            if vek_od is not None:
                podmienky.append(f"vek ≥ {vek_od}")
            if vek_do is not None:
                podmienky.append(f"vek ≤ {vek_do}")

            if podmienky:
                popis += " (" + ", ".join(podmienky) + ")"

            # chi-kvadrat: mutacia × pohlavie
            ct_pohlavie = pd.crosstab(df[mutacia], df["pohlavie"])
            if ct_pohlavie.shape[0] > 1 and ct_pohlavie.shape[1] > 1:
                observed_pohlavie = ct_pohlavie.values
                expected_pohlavie = chi2_contingency(ct_pohlavie)[3]
                degrees_of_freedom_pohlavie = (ct_pohlavie.shape[0] - 1) * (ct_pohlavie.shape[1] - 1)
                chi2, p, _ = chi_square_test(observed_pohlavie.flatten(),
                                             expected_pohlavie.flatten())
                vysledky.append(
                    f"{popis} × pohlavie: χ² = {chi2:.2f}, p = {p:.3f}, df = {degrees_of_freedom_pohlavie}")
            else:
                vysledky.append(f"{popis} × pohlavie: Nedostatok dát pre test.")

            # chi-kvadrat: mutacia × diagnoza
            if diagnoza and diagnoza.strip() != "" and diagnoza != "Všetky":
                diagnoza_maska = df["diagnoza MKCH-10"].apply(lambda x: diagnoza in str(x))
                ct_diagnoza = pd.crosstab(df[mutacia], diagnoza_maska)
                diagnoza_nazov = f"diagnóza '{diagnoza}'"
            else:
                ct_diagnoza = pd.crosstab(df[mutacia], pecen_maska)
                diagnoza_nazov = "pečeňová diagnóza ('K76.0', 'K75.9') "

            if ct_diagnoza.shape[0] > 1 and ct_diagnoza.shape[1] > 1:
                observed_diagnoza = ct_diagnoza.values
                expected_diagnoza = chi2_contingency(ct_diagnoza)[3]
                degrees_of_freedom_diagnoza = (ct_diagnoza.shape[0] - 1) * (ct_diagnoza.shape[1] - 1)
                chi2, p, _ = chi_square_test(observed_diagnoza.flatten(),
                                             expected_diagnoza.flatten())
                vysledky.append(
                    f"{popis} × {diagnoza_nazov}: χ² = {chi2:.2f}, p = {p:.3f}, df = {degrees_of_freedom_diagnoza}")
            else:
                vysledky.append(f"{popis} × {diagnoza_nazov}: Nedostatok dát pre test.")

        return ui.tags.div(*(ui.tags.p(v) for v in vysledky))

    @reactive.Effect
    def aktualizuj_vystupy():
        vek_od = input.vek_od()
        vek_do = input.vek_do()
        pohlavie = input.pohlavie()
        diagnoza = input.diagnoza()

        print(f"Rozmery df pred filtrovaním: {df.shape}")  # Debug
        filtrovane_data = priprav_data_pre_grafy(df, vek_od, vek_do, pohlavie, diagnoza)
        print(f"Rozmery filtrovaných dát: {filtrovane_data.shape}")  # Debug
        filtered_data.set(filtrovane_data)

    @output
    @render.ui
    def grafy_vystup():
        data = filtered_data.get()
        return generuj_grafy(data,
                             pohlavie=input.pohlavie(),
                             diagnoza=input.diagnoza(),
                             vek_od=input.vek_od(),
                             vek_do=input.vek_do())

    @output
    @render.ui
    def chi_kvadrat_vystup():
        data = filtered_data.get()
        return vykonaj_chi_kvadrat_testy(data, diagnoza=input.diagnoza(), vek_od=input.vek_od(), vek_do=input.vek_do())

    @output
    @render.ui
    def page_ui():
        if input.page() == "Úvod":
            return ui.TagList(
                ui.p(
                    "Vitajte v aplikácií, ktorá pomáha lekárom a výskumníkom analyzovať genetické dáta pacientov a ich diagnózy."),
                ui.p("Kľúčové funkcie:"),
                ui.p(
                    "   - Obsahuje očistený dataset : Opravené chyby a nepresnosti v údajoch (napr. chýbajúce hodnoty)."),
                ui.p(
                    "   - Analyzuje gény: Kontroluje, či rozloženie HFE génov je v norme (Hardy-Weinbergova rovnováha)."),
                ui.p(
                    "   - Určuje predispozície: Zisťuje, koľko pacientov má gény pre hemochromatózu (dedičnú chorobu)."),
                ui.p("   - Hľadá súvislosti: Skúma, či existuje spojitosť medzi génmi a chorobami, hlavne pečeňovými."),
                ui.p("   - Vytvára grafy: Zobrazuje dáta o génoch, veku, pohlaví a chorobách v prehľadných grafoch."),
                ui.p("   - Analyzuje diagnózy: Roztrieďuje diagnózy podľa MKCH-10 a sleduje, ako sa menili v čase."),
                ui.p("   - Overuje kódy: Kontroluje, či sú MKCH-10 kódy diagnóz správne a aktuálne.")
            )
        elif input.page() == "Hardy-Weinberg":
            return ui.TagList(
                ui.h2("Hardy-Weinbergova rovnováha"),
                ui.h3("Informácie o alelách a ich frekvencie"),
                ui.output_table("allele_info_table"),
                ui.h3("Pozorované počty genotypov v datasete"),
                ui.output_table("observed_values_table"),
                ui.h3("Očakávané počty genotypov podľa Hardy-Weinbergovej rovnováhy"),
                ui.output_table("expected_values_table"),
                ui.h3("Chi-kvadrát test"),
                ui.output_table("chi2_test_table"),
                ui.h3("Výsledky testu Hardy-Weinbergovej rovnováhy"),
                ui.output_table("hw_hypothesis_table"),
                ui.output_table("hw_results_table")
            )
        elif input.page() == "Genotypy a predispozície":
            results = analyze_genotype_distribution(df)
            total_patients = len(df)
            c282y_heterozygot_percent = results["genotype_percentages"].loc[
                results["genotype_percentages"]["Mutácia"] == "C282Y", "Heterozygot (%)"].iloc[0] if not results[
                "genotype_percentages"].empty and ("C282Y" in results["genotype_percentages"][
                "Mutácia"].values) else "N/A"
            s65c_homozygot_mutant_percent = results["genotype_percentages"].loc[
                results["genotype_percentages"]["Mutácia"] == "S65C", "Mutant (%)"].iloc[0] if not results[
                "genotype_percentages"].empty and ("S65C" in results["genotype_percentages"][
                "Mutácia"].values) else "N/A"
            compound_heterozygot_h63d_percent = results["risk_percentages"].loc[
                results["risk_percentages"]["Kategória"] == "C282Y/H63D zložený heterozygot", "Percento (%)"].iloc[
                0] if not results["risk_percentages"].empty and (
                    "C282Y/H63D zložený heterozygot" in results["risk_percentages"]["Kategória"].values) else "N/A"
            compound_heterozygot_c282y_s65c_percent = results["risk_percentages"].loc[
                results["risk_percentages"]["Kategória"] == "C282Y/S65C zložený heterozygot", "Percento (%)"].iloc[
                0] if not results["risk_percentages"].empty and (
                    "C282Y/S65C zložený heterozygot" in results["risk_percentages"]["Kategória"].values) else "N/A"
            h63d_s65c_compound_heterozygot_percent = "0.0"  # V dátach sa nevyskytuje

            return ui.TagList(
                ui.h2("Genotypy a predispozície k hemochromatóze"),
                ui.p("Analýza zastúpenia genotypov a predispozícií k hereditárnej hemochromatóze."),
                ui.h3("Percentuálne zastúpenie genotypov"),
                ui.output_table("genotype_distribution_table"),
                ui.h3("Predispozícia k hereditárnej hemochromatóze"),
                ui.output_table("hemochromatosis_risk_table"),
                ui.h3("Počet pacientov podľa predispozície"),
                ui.output_table("predisposition_summary_table")
            )
        elif input.page() == "Analýza diagnóz":
            vyskyt = analyzuj_diagnozy(df, "diagnoza MKCH-10", "validovany vysledok", mkch10_data)
            return ui.TagList(
                ui.h2("Analýza diagnóz podľa MKCH-10"),
                ui.output_ui("analyza_diagnoz_ui"),  # s podfarbenim ale html
                # ui.output_table("vyskyt_diagnoz_table") #bez podfarbenia ale v shiny style
            )
        elif input.page() == "MKCH-10":
            if not mkch10_data:
                return ui.p("Číselník MKCH-10 nebol načítaný.")

            names = sheet_names()

            return ui.TagList(
                ui.h2("MKCH-10 Číselník"),
                ui.div(
                    ui.div(
                        ui.input_text("mkch10_hladaj", "Hľadať v kóde alebo názve:"),
                        ui.div(
                            ui.input_action_button("mkch10_hladaj_button", "Hľadať", class_="action-button-sm"),
                            ui.input_action_button("mkch10_reset_search", "Zobraziť všetky", class_="action-button-sm"),
                            class_="horizontal-buttons"
                        ),
                        class_="vertical-layout align-items-start"
                    ),
                    ui.div(ui.output_text("mkch10_rows_count"), style="text-align: right;"),
                    class_="horizontal-layout justify-content-space-between align-items-start"
                ),
                ui.output_ui("mkch10_navigation"),
                ui.output_ui("mkch10_current_table")
            )

        elif input.page() == "Genotypy, demografia a diagnózy":
            return ui.TagList(
                ui.h2("Vzťah genotypov, demografie a diagnóz"),

                ui.input_numeric("vek_od", "Vek od:", value=0, min=0, max=120),
                ui.input_numeric("vek_do", "Vek do:", value=100, min=0, max=120),
                ui.input_select("pohlavie", "Pohlavie:", choices=["Všetky", "Muž", "Žena"], selected="Všetky"),
                ui.input_select("diagnoza", "Diagnóza (MKCH-10 kód):", choices=diagnozy, selected="Všetky"),

                ui.output_ui("grafy_vystup"),
                ui.output_ui("chi_kvadrat_vystup")
            )
        return ui.TagList(
            ui.h2("Očistený dataset"),
            ui.output_table("data_table")
        )

    @reactive.Effect
    @reactive.event(input.mkch10_hladaj_button)
    def _perform_search():
        hladany_vyraz = input.mkch10_hladaj()
        if mkch10_data:
            temp_filtered_data = {}
            sheet_names_list = list(mkch10_data.keys())
            for sheet_name, df in mkch10_data.items():
                if sheet_name == nazvy_harkov_mkch10[0]:
                    continue

                kod_col = 'Kód diagnózy' if 'Kód diagnózy' in df.columns else 'Kód' if 'Kód' in df.columns else None
                nazov_col = 'Nazov' if 'Nazov' in df.columns else 'Názov' if 'Názov' in df.columns else None

                if kod_col and nazov_col:
                    filtered_df = df[
                        df[kod_col].str.contains(hladany_vyraz, case=False, na=False) |
                        df[nazov_col].str.contains(hladany_vyraz, case=False, na=False)
                        ].fillna('')
                    if not filtered_df.empty:
                        temp_filtered_data[sheet_name] = filtered_df
                else:
                    logging.warning(f"V hárku '{sheet_name}' sa nenašli stĺpce pre kód alebo názov.")

            filtered_mkch10_data.set(temp_filtered_data)
            search_performed.set(True)
            current_sheet_index.set(0)
        else:
            filtered_mkch10_data.set(None)
            search_performed.set(False)

    @output
    @render.ui
    def mkch10_current_table():
        current_index = current_sheet_index()
        rendered_table = ui.div("Žiadne výsledky vyhľadávania.")

        if search_performed.get() and filtered_mkch10_data.get():
            filtered_data = list(filtered_mkch10_data().values())
            filtered_sheet_names = list(filtered_mkch10_data().keys())

            if filtered_data:
                if current_index < len(filtered_data):
                    data_to_render = filtered_data[current_index]
                    current_sheet_name = filtered_sheet_names[current_index]

                    def format_row(row):
                        podfarbit = False
                        kod_col_name = 'Kód' if 'Kód' in data_to_render.columns else 'Kód diagnózy' if 'Kód diagnózy' in data_to_render.columns else None
                        if kod_col_name and '-' in str(row[kod_col_name]):
                            podfarbit = True
                        cells = [tags.td("" if pd.isna(value) else str(value)) for value in row.values]
                        return tags.tr(style="background-color: #F0F0F0;" if podfarbit else "", *cells)

                    table_header = tags.thead(
                        tags.tr(*[tags.th(col, style="background-color: #F0F8FF;") for col in data_to_render.columns]))
                    table_body = tags.tbody(*[format_row(row) for index, row in data_to_render.iterrows()])
                    rendered_table = ui.TagList(
                        ui.div(f"Výsledky z hárku: {current_sheet_name}",
                               style="font-size: 1.1em; font-weight: bold; margin-bottom: 5px;"),
                        tags.table(table_header, table_body, class_="dataframe")
                    )
            else:
                rendered_table = ui.div("Žiadne výsledky pre zadaný výraz.")

        elif not search_performed.get():
            current_sheet_name = sheet_names()[current_index]
            data_to_render = mkch10_data.get(current_sheet_name)
            if data_to_render is not None:
                def format_row(row):
                    podfarbit = False
                    kod_col_name = 'Kód' if 'Kód' in data_to_render.columns else 'Kód diagnózy' if 'Kód diagnózy' in data_to_render.columns else None
                    if kod_col_name and '-' in str(row[kod_col_name]):
                        podfarbit = True
                    cells = [tags.td("" if pd.isna(value) else str(value)) for value in row.values]
                    return tags.tr(style="background-color: #F0F0F0;" if podfarbit else "", *cells)

                table_header = tags.thead(
                    tags.tr(*[tags.th(col, style="background-color: #F0F8FF;") for col in data_to_render.columns]))
                table_body = tags.tbody(*[format_row(row) for index, row in data_to_render.iterrows()])
                rendered_table = tags.table(table_header, table_body, class_="dataframe")
            else:
                rendered_table = ui.div("Dáta pre tento hárok neboli načítané.")

        elif search_performed.get() and not filtered_mkch10_data.get():
            rendered_table = ui.div("Žiadne výsledky pre zadaný výraz.")

        return rendered_table

    @output
    @render.ui
    def mkch10_navigation():
        names = sheet_names()
        nav_buttons = []
        import random
        if search_performed.get() and filtered_mkch10_data.get():
            filtered_sheet_names = list(filtered_mkch10_data().keys())
            for i, name in enumerate(filtered_sheet_names):
                nav_buttons.append(ui.input_action_button(f"goto_sheet_filtered_{i}", name, class_="action-button-nav",
                                                          style="margin-right: 5px;",
                                                          key=f"filtered_nav_{i}_{random.random()}"))
        else:
            for i, name in enumerate(names):
                nav_buttons.append(ui.input_action_button(f"goto_sheet_{i}", name, class_="action-button-nav",
                                                          style="margin-right: 5px;", key=f"nav_{i}_{random.random()}"))
        return ui.div(*nav_buttons, style="margin-bottom: 10px; overflow-x: auto; white-space: nowrap;")

    previous_clicks = {}

    @reactive.Effect
    def _update_sheet_index_from_nav():
        names = sheet_names()
        if search_performed.get() and filtered_mkch10_data.get():
            filtered_sheet_names = list(filtered_mkch10_data().keys())
            for i in range(len(filtered_sheet_names)):
                btn_id = f"goto_sheet_filtered_{i}"
                current_clicks = input[btn_id]()
                previous = previous_clicks.get(btn_id, 0)
                if current_clicks > previous:
                    print(f"Kliknuté (value) na filtrované tlačidlo s indexom: {i}")
                    current_sheet_index.set(i)
                previous_clicks[btn_id] = current_clicks
        else:
            for i in range(len(names)):
                btn_id = f"goto_sheet_{i}"
                current_clicks = input[btn_id]()
                previous = previous_clicks.get(btn_id, 0)
                if current_clicks > previous:
                    print(f"Kliknuté (value) na tlačidlo s indexom: {i}")
                    current_sheet_index.set(i)
                previous_clicks[btn_id] = current_clicks

    @reactive.Effect
    @reactive.event(input.mkch10_reset_search)
    def _reset_search():
        search_performed.set(False)
        current_sheet_index.set(1)

    @output
    @render.text
    def mkch10_rows_count():
        if search_performed.get() and filtered_mkch10_data.get():
            filtered_data = list(filtered_mkch10_data().values())
            if filtered_data and current_sheet_index() < len(filtered_data):
                return f"Počet nájdených riadkov: {len(filtered_data[current_sheet_index()])}"
            else:
                return "Žiadne výsledky."
        elif mkch10_data and sheet_names():
            current_sheet_name = sheet_names()[current_sheet_index()]
            data_to_render = mkch10_data.get(current_sheet_name)
            if data_to_render is not None:
                return f"Počet riadkov v aktuálnom hárku: {len(data_to_render)}"
            else:
                return ""
        else:
            return ""

    @output
    @render.table
    def data_table():
        return df

    def generate_hw_table(df, value_extractor):
        results = []
        for column in ["HFE C187G (H63D) [HFE]", "HFE A193T (S65C) [HFE]", "HFE G845A (C282Y) [HFE]"]:
            hw_result = check_hardy_weinberg(df, column)
            if hw_result is not None:
                results.append(value_extractor(hw_result, column))
            else:
                results.append(value_extractor({
                    "Mutácia": column,
                    "reason_na": "Chyba pri výpočte (pravdepodobne nedostatok dát)",
                    "allele_info": None,
                    "observed": None,
                    "expected": None,
                    "chi2_results": None
                }, column))
        final_df = pd.DataFrame(results)
        if 'Dôvod N/A' in final_df.columns:
            non_default_reasons = final_df[final_df['Dôvod N/A'].notna() & (final_df['Dôvod N/A'] != 'Dostatočné dáta')]
            if non_default_reasons.empty:
                final_df = final_df.drop(columns=['Dôvod N/A'])
        return final_df

    @output
    @render.table
    def allele_info_table():
     return generate_hw_table(df, lambda r, c: {
         "Mutácia": r.get("Mutácia", c),
         "Celkový počet alel": r.get("allele_info", {}).get("total_alleles") if r and r.get("allele_info") is not None else "N/A",
         "Počet normálnych alel": r.get("allele_info", {}).get("normal_alleles") if r and r.get("allele_info") is not None else "N/A",
         "Počet mutovaných alel": r.get("allele_info", {}).get("mutant_alleles") if r and r.get("allele_info") is not None else "N/A",
         "Frekvencia normálnej alely (p)": f"{r.get('allele_info', {}).get('allele_p', 0):.5f}" if r and r.get("allele_info") is not None else "N/A",
         "Frekvencia mutovanej alely (q)": f"{r.get('allele_info', {}).get('allele_q', 0):.5f}" if r and r.get("allele_info") is not None else "N/A",
         "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
     })

    @output
    @render.table
    def observed_values_table():
     return generate_hw_table(df, lambda r, c: {
         "Mutácia": r.get("Mutácia", c),
         "Normal": r.get("observed", {}).get("normal") if r and r.get("observed") is not None else "N/A",
         "Heterozygot": r.get("observed", {}).get("heterozygot") if r and r.get("observed") is not None else "N/A",
         "Mutant": r.get("observed", {}).get("mutant") if r and r.get("observed") is not None else "N/A",
         "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
     })

    @output
    @render.table
    def expected_values_table():
     return generate_hw_table(df, lambda r, c: {
         "Mutácia": r.get("Mutácia", c),
         "Normal": round(r.get("expected", {}).get("normal"), 2) if r and r.get("expected") is not None and r.get("expected").get("normal") is not None else "N/A",
         "Heterozygot": round(r.get("expected", {}).get("heterozygot"), 2) if r and r.get("expected") is not None and r.get("expected").get("heterozygot") is not None else "N/A",
         "Mutant": round(r.get("expected", {}).get("mutant"), 2) if r and r.get("expected") is not None and r.get("expected").get("mutant") is not None else "N/A",
         "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
     })

    @output
    @render.table
    def chi2_test_table():
     return generate_hw_table(df, lambda r, c: {
         "Mutácia": r.get("Mutácia", c),
         "Chi-kvadrát test": round(r.get("chi2_results", {}).get("chi2"), 5) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("chi2") is not None else "N/A",
         "Stupne voľnosti": r.get("chi2_results", {}).get("df") if r and r.get("chi2_results") is not None else "N/A",
         "p-hodnota": format_p_value(r.get("chi2_results", {}).get("p_value")) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None else "N/A",
         "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
     })

    @output
    @render.table
    def hw_hypothesis_table():
     hypothesis_df = pd.DataFrame([
         {
             "Hypotézy": "H0: Distribúcia HFE genotypov pre mutáciu je v Hardy-Weinbergovej rovnováhe (očakávané frekvencie genotypov zodpovedajú pozorovaným)",
         },
         {
             "Hypotézy": "H1: Distribúcia HFE genotypov pre mutáciu nie je v Hardy-Weinbergovej rovnováhe (pozorované frekvencie genotypov sa významne líšia od očakávaných)",
         }
     ])
     return hypothesis_df

    @output
    @render.table
    def hw_results_table():
     return generate_hw_table(df, lambda r, c: {
         "Mutácia": r.get("Mutácia", c),
         "Chi-kvadrát test": round(r.get("chi2_results", {}).get("chi2"), 5) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("chi2") is not None else "N/A",
         "p-hodnota": format_p_value(r.get("chi2_results", {}).get("p_value")) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None else "N/A",
         "Výsledok": (
             "Na základe chi-kvadrát testu sa nulová hypotéza H0 o Hardy-Weinbergovej rovnováhe nezamieta."
             if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None and r["chi2_results"]["p_value"] > 0.05
             else (
                 "Na základe chi-kvadrát testu sa nulová hypotéza H0 o Hardy-Weinbergovej rovnováhe zamieta a prijíma sa R1 hypotéza."
                 if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None
                 else "N/A"
             )
         ),
         "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
     })


    @output
    @render.ui
    def analyza_diagnoz_ui():
        analyza = analyzuj_diagnozy(df, "diagnoza MKCH-10", "validovany vysledok", mkch10_data)
        vyskyt_diagnoz_df = analyza['vyskyt_diagnoz']
        celkovy_pocet = analyza['celkovy_pocet']
        chybne_kody_df = analyza['chybne_kody']

        def get_row_style(row):
            return 'background-color: lightcoral;' if bool(row.get('Je_chybny', False)) else ''

        table_html = "<table>"
        table_html += "<thead><tr><th>" + "</th><th>".join(vyskyt_diagnoz_df.columns) + "</th></tr></thead>"
        table_html += "<tbody>"
        for index, row in vyskyt_diagnoz_df.iterrows():
            style = get_row_style(row)
            table_html += f"<tr style='{style}'>"
            table_html += "<td>" + "</td><td>".join(map(str, row.values)) + "</td>"
            table_html += "</tr>"
        table_html += "</tbody></table>"

        return ui.TagList(
            ui.p(f"Celkový počet záznamov: {celkovy_pocet}"),
            ui.HTML(table_html),
            ui.h4("Chybné kódy"),
            ui.p(", ".join(
                chybne_kody_df["diagnoza MKCH-10"].tolist()) if not chybne_kody_df.empty else "Žiadne chybné kódy")
        )

    @output
    @render.text
    def chybne_kody_text():
        analyza = analyzuj_diagnozy(df, "diagnoza MKCH-10", "validovany vysledok", mkch10_data)
        chybne_kody = analyza['chybne_kody']
        return f"Chybné kódy: {', '.join(chybne_kody['diagnoza MKCH-10'])}" if not chybne_kody.empty else "Žiadne chybné kódy"

    @output
    @render.table
    def genotype_distribution_table():
        return analyze_genotype_distribution(df)["genotype_percentages"]

    @output
    @render.table
    def hemochromatosis_risk_table():
        return analyze_genotype_distribution(df)["risk_percentages"]

    @output
    @render.table
    def predisposition_summary_table():
     predisposition_df = analyze_genotype_distribution(df)["predisposition_table"]
     predisposition_df.loc[predisposition_df["Skupina"] == "Celkový počet pacientov", "Percento (%)"] = ""
     return predisposition_df

    @output
    @render.table
    def vyskyt_diagnoz_table():
     analyza = analyzuj_diagnozy(df, "diagnoza MKCH-10", "validovany vysledok", mkch10_data)
     return analyza['vyskyt_diagnoz']


app = App(app_ui, server)
