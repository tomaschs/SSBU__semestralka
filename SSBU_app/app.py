from shiny import App, render, ui, reactive
from shiny.ui import tags
import pandas as pd
from shared import df, mkch10_file_path, nazvy_harkov_mkch10
from app_ui import app_ui
from utils import check_hardy_weinberg, analyze_genotype_distribution, analyzuj_diagnozy, prirad_kapitolu_mkch10, nacitaj_mkch10_ciselnik

# Načítame číselník pri spustení aplikácie
mkch10_data = nacitaj_mkch10_ciselnik(str(mkch10_file_path), nazvy_harkov_mkch10)

def server(input, output, session):
    sheet_names = reactive.Value(list(mkch10_data.keys()) if mkch10_data else [])
    current_sheet_index = reactive.Value(0)
    filtered_mkch10_data = reactive.Value(None) # Bude obsahovať filtrovaný DataFrame
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
                ui.h3("Percentuálne zastúpenie špecifických genotypov"),
                ui.p(f"Percentuálne zastúpenie C282Y heterozygotov: {c282y_heterozygot_percent}%"),
                ui.p(f"Percentuálne zastúpenie C282Y/H63D zložených heterozygotov: {compound_heterozygot_h63d_percent}%"),
                ui.p(f"Percentuálne zastúpenie C282Y/S65C zložených heterozygotov: {compound_heterozygot_c282y_s65c_percent}% (v dátach 0.0%)"),
                ui.p(f"Percentuálne zastúpenie H63D/S65C zložených heterozygotov: {h63d_s65c_compound_heterozygot_percent}% (v dátach 0.0%)"),
                ui.p(f"Percentuálne zastúpenie S65C homozygotných mutantov: {s65c_homozygot_mutant_percent}%"),
                ui.h3("Percentuálne zastúpenie genotypov"),
                ui.output_table("genotype_distribution_table"),
                ui.h3("Predispozícia k hereditárnej hemochromatóze"),
                ui.output_table("hemochromatosis_risk_table"),
                ui.h3("Počet pacientov podľa predispozície"),
                ui.output_table("predisposition_summary_table")
            )
        elif input.page() == "Analýza diagnóz":
            vyskyt = analyzuj_diagnozy(df, "diagnoza MKCH-10", "validovany vysledok")
            return ui.TagList(
                ui.h2("Analýza diagnóz podľa MKCH-10"),
                ui.output_table("vyskyt_diagnoz_table"),
            )
        elif input.page() == "MKCH-10":
            if not mkch10_data:
                return ui.p("Číselník MKCH-10 nebol načítaný.")

            names = sheet_names()
            current_index = current_sheet_index()

            return ui.TagList(
                ui.h2("MKCH-10 Číselník"),
                ui.input_text("mkch10_hladaj", "Hľadať v kóde alebo názve:"),
                ui.input_action_button("mkch10_hladaj_button", "Hľadať"),
                ui.div(names[current_index], style="font-size: 1.5em; font-weight: bold; margin-bottom: 10px;"),
                ui.output_ui("mkch10_current_table"),
                ui.div(
                    ui.input_action_button("prev_sheet", "Predchádzajúci"),
                    ui.span(" "),
                    ui.input_action_button("next_sheet", "Nasledujúci")
                )
            )
        else:
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
            for i, sheet_name in enumerate(sheet_names_list):
                if i == 0:
                    continue

                df = mkch10_data[sheet_name]
                kod_col = 'Kód diagnózy' if 'Kód diagnózy' in df.columns else 'Kód' if 'Kód' in df.columns else None
                nazov_col = 'Nazov' if 'Nazov' in df.columns else 'Názov' if 'Názov' in df.columns else None

                if kod_col and nazov_col:
                    filtered_df = df[
                        df[kod_col].str.contains(hladany_vyraz, case=False, na=False) |
                        df[nazov_col].str.contains(hladany_vyraz, case=False, na=False)
                    ].fillna('')
                    temp_filtered_data[sheet_name] = filtered_df
                else:
                    temp_filtered_data[sheet_name] = pd.DataFrame() # Ak stĺpce nenájdené, prázdny DataFrame
            filtered_mkch10_data.set(temp_filtered_data)
            search_performed.set(True)
        else:
            filtered_mkch10_data.set(None)
            search_performed.set(False)

    @output
    @render.ui
    def mkch10_current_table():
        current_index = current_sheet_index()
        current_sheet_name = sheet_names()[current_index]
        data_to_render = mkch10_data.get(current_sheet_name)

        if search_performed.get() and filtered_mkch10_data.get() is not None:
            data_to_render = filtered_mkch10_data.get().get(current_sheet_name, pd.DataFrame())

        if data_to_render is not None:
            def format_row(row):
                podfarbit = False
                kod_col_name = 'Kód' if 'Kód' in data_to_render.columns else 'Kód diagnózy' if 'Kód diagnózy' in data_to_render.columns else None
                if kod_col_name and '-' in str(row[kod_col_name]):
                    podfarbit = True
                cells = [tags.td(value) for value in row.values.astype(str)]
                if podfarbit:
                    return tags.tr(style="background-color: #F0F0F0;", *cells)
                else:
                    return tags.tr(*cells)

            table_header = tags.thead(tags.tr(*[tags.th(col, style="background-color: #F0F8FF;") for col in data_to_render.columns]))
            table_body = tags.tbody(*[format_row(row) for index, row in data_to_render.iterrows()])
            return tags.table(table_header, table_body, class_="dataframe")
        else:
            return ui.div("Dáta MKCH-10 neboli načítané.")

    @reactive.Effect
    @reactive.event(input.prev_sheet)
    def _prev_sheet():
        search_performed.set(False) # Reset vyhľadávanie pri zmene hárku
        current_sheet_index.set(max(0, current_sheet_index() - 1))

    @reactive.Effect
    @reactive.event(input.next_sheet)
    def _next_sheet():
        search_performed.set(False) # Reset vyhľadávanie pri zmene hárku
        current_sheet_index.set(min(len(sheet_names()) - 1, current_sheet_index() + 1))

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
            "Frekvencia normálnej alely (p)": f"{r.get('allele_info', {}).get('allele_p', 0):.4f}" if r and r.get("allele_info") is not None else "N/A",
            "Frekvencia mutovanej alely (q)": f"{r.get('allele_info', {}).get('allele_q', 0):.4f}" if r and r.get("allele_info") is not None else "N/A",
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
            "Chi-kvadrát test": round(r.get("chi2_results", {}).get("chi2"), 4) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("chi2") is not None else "N/A",
            "Stupne voľnosti": r.get("chi2_results", {}).get("df") if r and r.get("chi2_results") is not None else "N/A",
            "p-hodnota": round(r.get("chi2_results", {}).get("p_value"), 4) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None else "N/A",
            "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
        })

    @output
    @render.table
    def hw_results_table():
        return generate_hw_table(df, lambda r, c: {
            "Mutácia": r.get("Mutácia", c),
            "Chi-kvadrát test": round(r.get("chi2_results", {}).get("chi2"), 4) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("chi2") is not None else "N/A",
            "p-hodnota": round(r.get("chi2_results", {}).get("p_value"), 4) if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None else "N/A",
            "Výsledok": (
                "Populácia je v Hardy-Weinbergovej rovnováhe (očakávané frekvencie genotypov zodpovedajú pozorovaným)."
                if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None and r["chi2_results"]["p_value"] > 0.05
                else (
                    "Populácia nie je v Hardy-Weinbergovej rovnováhe (pozorované frekvencie genotypov sa významne líšia od očakávaných), čo môže naznačovať vplyv evolučných síl alebo iných faktorov."
                    if r and r.get("chi2_results") is not None and r.get("chi2_results").get("p_value") is not None
                    else "N/A"
                )
            ),
            "Dôvod N/A": r.get("reason_na", "Dostatočné dáta") if r is not None else "Chyba pri výpočte (pravdepodobne nedostatok dát)"
        })

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
        return analyze_genotype_distribution(df)["predisposition_table"]

    @output
    @render.table
    def vyskyt_diagnoz_table():
        return analyzuj_diagnozy(df, "diagnoza MKCH-10", "validovany vysledok")

app = App(app_ui, server)