import logging
from shiny import App, render, ui, reactive
from shiny.ui import tags
import pandas as pd
from shared import df, mkch10_file_path, nazvy_harkov_mkch10
from app_ui import app_ui
from utils import check_hardy_weinberg, analyze_genotype_distribution, analyzuj_diagnozy, prirad_kapitolu_mkch10, nacitaj_mkch10_ciselnik, format_p_value

mkch10_data = nacitaj_mkch10_ciselnik(str(mkch10_file_path), nazvy_harkov_mkch10)

def server(input, output, session):
     sheet_names = reactive.Value(list(mkch10_data.keys()) if mkch10_data else [])
     current_sheet_index = reactive.Value(0)
     filtered_mkch10_data = reactive.Value(None)
     search_performed = reactive.Value(False)

     @output
     @render.ui
     def page_ui():
         if input.page() == "Úvod":
             return ui.TagList(
                 ui.h2("Vitajte v SSBU aplikacii"),
                 ui.p("Úvodný obsah aplikácie...")
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
                 ui.output_table("hw_results_table")
             )
         elif input.page() == "Genotypy a predispozície":
             results = analyze_genotype_distribution(df)
             total_patients = len(df)
             c282y_heterozygot_percent = results["genotype_percentages"].loc[results["genotype_percentages"]["Mutácia"] == "C282Y", "Heterozygot (%)"].iloc[0] if not results["genotype_percentages"].empty and ("C282Y" in results["genotype_percentages"]["Mutácia"].values) else "N/A"
             s65c_homozygot_mutant_percent = results["genotype_percentages"].loc[results["genotype_percentages"]["Mutácia"] == "S65C", "Mutant (%)"].iloc[0] if not results["genotype_percentages"].empty and ("S65C" in results["genotype_percentages"]["Mutácia"].values) else "N/A"
             compound_heterozygot_h63d_percent = results["risk_percentages"].loc[results["risk_percentages"]["Kategória"] == "C282Y/H63D zložený heterozygot", "Percento (%)"].iloc[0] if not results["risk_percentages"].empty and ("C282Y/H63D zložený heterozygot" in results["risk_percentages"]["Kategória"].values) else "N/A"
             compound_heterozygot_c282y_s65c_percent = results["risk_percentages"].loc[results["risk_percentages"]["Kategória"] == "C282Y/S65C zložený heterozygot", "Percento (%)"].iloc[0] if not results["risk_percentages"].empty and ("C282Y/S65C zložený heterozygot" in results["risk_percentages"]["Kategória"].values) else "N/A"
             h63d_s65c_compound_heterozygot_percent = "0.0" # V dátach sa nevyskytuje

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
                 ui.output_ui("analyza_diagnoz_ui"), #s podfarbenim ale html
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
                        ui.div( # Nový div pre tlačidlá
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
             for sheet_name, df in mkch10_data.items(): # Iterujeme priamo cez názov hárku a DataFrame
                 if sheet_name == nazvy_harkov_mkch10[0]: # Preskočí prvý hárok s vysvetlivkami
                     continue

                 kod_col = 'Kód diagnózy' if 'Kód diagnózy' in df.columns else 'Kód' if 'Kód' in df.columns else None
                 nazov_col = 'Nazov' if 'Nazov' in df.columns else 'Názov' if 'Názov' in df.columns else None

                 if kod_col and nazov_col:
                     filtered_df = df[
                         df[kod_col].str.contains(hladany_vyraz, case=False, na=False) |
                         df[nazov_col].str.contains(hladany_vyraz, case=False, na=False)
                     ].fillna('')
                     if not filtered_df.empty: # Ukladáme len neprázdne výsledky
                         temp_filtered_data[sheet_name] = filtered_df
                 else:
                     logging.warning(f"V hárku '{sheet_name}' sa nenašli stĺpce pre kód alebo názov.")

             filtered_mkch10_data.set(temp_filtered_data)
             search_performed.set(True)
             current_sheet_index.set(0) # Po vyhľadávaní nastavíme na prvý výsledok
         else:
             filtered_mkch10_data.set(None)
             search_performed.set(False)

     @output
     @render.ui
     def mkch10_current_table():
            current_index = current_sheet_index()
            rendered_table = ui.div("Žiadne výsledky vyhľadávania.")  # Predvolený text

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

                        table_header = tags.thead(tags.tr(*[tags.th(col, style="background-color: #F0F8FF;") for col in data_to_render.columns]))
                        table_body = tags.tbody(*[format_row(row) for index, row in data_to_render.iterrows()])
                        rendered_table = ui.TagList(
                            ui.div(f"Výsledky z hárku: {current_sheet_name}", style="font-size: 1.1em; font-weight: bold; margin-bottom: 5px;"),
                            tags.table(table_header, table_body, class_="dataframe")
                        )
                else:
                    rendered_table = ui.div("Žiadne výsledky pre zadaný výraz.")
            
            elif not search_performed.get():  # Ak nebolo vyhľadávané, zobrazíme aktuálny hárok z pôvodných dát
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

                    table_header = tags.thead(tags.tr(*[tags.th(col, style="background-color: #F0F8FF;") for col in data_to_render.columns]))
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
                nav_buttons.append(ui.input_action_button(f"goto_sheet_filtered_{i}", name, class_="action-button-nav", style="margin-right: 5px;", key=f"filtered_nav_{i}_{random.random()}"))
        else:
            for i, name in enumerate(names):
                nav_buttons.append(ui.input_action_button(f"goto_sheet_{i}", name, class_="action-button-nav", style="margin-right: 5px;", key=f"nav_{i}_{random.random()}"))
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
         current_sheet_index.set(1) # Návrat na prvý dátový hárok

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
     def hw_results_table():
         return generate_hw_table(df, lambda r, c: {
             "Mutácia": r.get("Mutácia", c),
             "Chi-kvadrát test": round(r.get("chi2_results", {}).get("chi2"), 5) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("chi2") is not None else "N/A",
             "p-hodnota": format_p_value(r.get("chi2_results", {}).get("p_value")) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None else "N/A",
             "Výsledok": (
                 "Populácia je v Hardy-Weinbergovej rovnováhe (očakávané frekvencie genotypov zodpovedajú pozorovaným)."
                 if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None and r["chi2_results"]["p_value"] > 0.05
                 else (
                     "Populácia nie je v Hardy-Weinbergovej rovnováhe (pozorované frekvencie genotypov sa významne líšia od očakávaných)."
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

        # Funkcia na generovanie štýlu riadku
        def get_row_style(row):
         return 'background-color: lightcoral;' if bool(row.get('Je_chybny', False)) else ''


        # Generovanie HTML tabuľky so štýlmi riadkov
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
            ui.p(", ".join(chybne_kody_df["diagnoza MKCH-10"].tolist()) if not chybne_kody_df.empty else "Žiadne chybné kódy")
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
