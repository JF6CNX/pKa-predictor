import flet as ft
from flet_webview import WebView
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import base64
import pubchempy as pcp
import re
import json
import os
import requests
from bs4 import BeautifulSoup

# --- 定数定義 ---
USER_DATA_FILE = "user_data.json"
PKa_RULEBOOK = {
    "sulfonic_acid": {"name": "スルホン酸", "pka": -3.0, "smarts": "S(=O)(=O)O", "type": "acid"},
    "carboxylic_acid": {"name": "カルボン酸", "pka": 4.8, "smarts": "[CX3](=O)[OX2H1]", "type": "acid"},
    "phenol": {"name": "フェノール", "pka": 10.0, "smarts": "c1ccccc1O", "type": "acid"},
    "alcohol": {"name": "アルコール", "pka": 16.0, "smarts": "[#6][OX2H]", "type": "acid"},
    "amine": {"name": "アミン", "pka": 10.6, "smarts": "[#7;!H0;!$(*C=O)]", "type": "base"},
    "pyridine": {"name": "ピリジン", "pka": 5.2, "smarts": "n1ccccc1", "type": "base"},
    "amide": {"name": "アミド", "pka": -0.5, "smarts": "[#7]C=O", "type": "base"},
}

# --- 起動時にSMARTSパターンをコンパイル ---
COMPILED_PATTERNS = {
    key: Chem.MolFromSmarts(rule["smarts"]) for key, rule in PKa_RULEBOOK.items()
}

# --- ユーザーデータを読み書きする関数 ---
def load_user_data():
    if os.path.exists(USER_DATA_FILE):
        with open(USER_DATA_FILE, 'r', encoding='utf-8') as f:
            try: return json.load(f)
            except json.JSONDecodeError: return {}
    return {}

def save_user_data(data):
    with open(USER_DATA_FILE, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=4, ensure_ascii=False)

def main(page: ft.Page):
    page.title = "pKa予測アプリ"
    page.window_width = 800
    page.window_height = 900
    page.theme_mode = ft.ThemeMode.LIGHT
    
    current_png_data = None
    current_predicted_pka = None

    # --- ご提供いただいた新しい物理化学データ取得関数 ---
    def get_physicochemical_data(compound):
        user_data = load_user_data()
        smiles = compound.smiles
        if smiles in user_data and all(k in user_data[smiles] for k in ["boiling_point", "melting_point", "density"]):
            return user_data[smiles]

        cid = compound.cid
        if not cid:
            return {"boiling_point": "データなし", "melting_point": "データなし", "density": "データなし"}

        def convert_to_celsius(value):
            """°F→°C変換（数値部分のみ）"""
            try:
                m = re.search(r"([0-9.-]+)", value)
                if not m:
                    return value
                num = float(m.group(1))
                if "°F" in value or "F" in value:
                    num = (num - 32) * 5 / 9
                return f"{num:.1f} °C"
            except:
                return value

        data_dict = {"boiling_point": "データなし", "melting_point": "データなし", "density": "データなし"}

        try:
            # --- PubChem APIから取得 ---
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/BoilingPoint,MeltingPoint,Density/JSON"
            r = requests.get(url, timeout=5)
            if r.status_code == 200:
                props = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
                if "BoilingPoint" in props:
                    data_dict["boiling_point"] = f"{props['BoilingPoint']} °C (API)"
                if "MeltingPoint" in props:
                    data_dict["melting_point"] = f"{props['MeltingPoint']} °C (API)"
                if "Density" in props:
                    data_dict["density"] = f"{props['Density']} g/cm³ (API)"
        except:
            pass

        # --- PUG View fallback ---
        if any(v == "データなし" for v in data_dict.values()):
            try:
                view_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
                r = requests.get(view_url, timeout=5)
                if r.status_code == 200:
                    record = r.json().get("Record", {})
                    sections = record.get("Section", [])
                    def walk(sec):
                        heading = sec.get("TOCHeading", "").lower()
                        infos = sec.get("Information", [])
                        for inf in infos:
                            val = inf.get("Value", {}).get("StringWithMarkup", [{}])[0].get("String", "")
                            if "boiling point" in heading and data_dict["boiling_point"] == "データなし":
                                data_dict["boiling_point"] = f"{convert_to_celsius(val)} (View)"
                            elif "melting point" in heading and data_dict["melting_point"] == "データなし":
                                data_dict["melting_point"] = f"{convert_to_celsius(val)} (View)"
                            elif "density" in heading and data_dict["density"] == "データなし":
                                data_dict["density"] = f"{val} (View)"
                        for s in sec.get("Section", []):
                            walk(s)
                    for s in sections:
                        walk(s)
            except:
                pass

        # --- Webスクレイピング fallback ---
        if any(v == "データなし" for v in data_dict.values()):
            try:
                html = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}", timeout=5).text
                soup = BeautifulSoup(html, "html.parser")
                text = soup.get_text("\n")
                for line in text.splitlines():
                    if "Boiling Point" in line and data_dict["boiling_point"] == "データなし":
                        data_dict["boiling_point"] = f"{convert_to_celsius(line.strip())} (Web)"
                    if "Melting Point" in line and data_dict["melting_point"] == "データなし":
                        data_dict["melting_point"] = f"{convert_to_celsius(line.strip())} (Web)"
                    if "Density" in line and data_dict["density"] == "データなし":
                        data_dict["density"] = f"{line.strip()} (Web)"
            except:
                pass

        user_data.setdefault(smiles, {}).update(data_dict)
        save_user_data(user_data)
        return data_dict


    # --- ここから下は既存のコード（一部修正） ---
    def save_file_result(e: ft.FilePickerResultEvent):
        if e.path and current_png_data:
            with open(e.path, "wb") as f: f.write(current_png_data)
            page.snack_bar = ft.SnackBar(ft.Text(f"画像を {e.path} に保存しました。"), duration=3000)
            page.snack_bar.open = True; page.update()
    save_dialog = ft.FilePicker(on_result=save_file_result)
    
    def open_context_menu(e: ft.TapEvent):
        if current_png_data and mol_image.visible:
            context_menu.left = e.global_x; context_menu.top = e.global_y
            context_menu.visible = True; page.update()

    def close_context_menu(e):
        if context_menu.visible:
            context_menu.visible = False; page.update()

    def save_image_click(e):
        close_context_menu(e)
        if current_png_data:
            save_dialog.save_file(dialog_title="画像を保存", file_name=f"{(id_textbox.value or 'molecule').replace(' ', '_')}.png", allowed_extensions=["png"])
        else:
            page.snack_bar = ft.SnackBar(ft.Text("保存する画像がありません。"), duration=3000, bgcolor="orange_accent_700")
            page.snack_bar.open = True; page.update()
    
    context_menu = ft.Card(content=ft.Column([ft.ListTile(title=ft.Text("画像を保存", color="black"), on_click=save_image_click, autofocus=True, hover_color="grey_200", width=200)], spacing=0),
        elevation=4, visible=False, color="white", shape=ft.RoundedRectangleBorder(radius=ft.border_radius.all(8)))
    
    def on_candidate_click(e):
        page.dialog.open = False
        run_prediction(e.control.data)
    
    candidate_dialog = ft.AlertDialog(modal=True, title=ft.Text("複数の候補が見つかりました"), content=ft.Column(spacing=5, height=200, scroll=ft.ScrollMode.ADAPTIVE))
    page.overlay.extend([save_dialog, context_menu, candidate_dialog]); page.on_click = close_context_menu

    def find_cas_number(synonyms):
        cas_pattern = re.compile(r'^\d{2,7}-\d{2}-\d$')
        for syn in synonyms:
            if cas_pattern.match(syn): return syn
        return "N/A"
    
    def resolve_identifier_to_smiles(identifier, id_type):
        try:
            namespace = 'name' if id_type == 'cas' else id_type
            compounds = pcp.get_compounds(identifier, namespace)
            if compounds: return compounds, None
        except: pass
        if id_type == 'name' or id_type == 'cas':
            try:
                url = f"https://cactus.nci.nih.gov/chemical/structure/{identifier}/smiles"
                response = requests.get(url, timeout=5)
                if response.status_code == 200 and response.text.strip():
                    compounds = pcp.get_compounds(response.text.strip(), 'smiles')
                    if compounds: return compounds, None
            except: pass
        return [], f"'{identifier}' は見つかりませんでした。"

    def predict_pka(smiles: str):
        nonlocal current_png_data, current_predicted_pka
        current_predicted_pka, current_png_data = None, None
        try: mol = Chem.MolFromSmiles(smiles)
        except: return "無効なSMILESです。", None, None, None
        if not mol: return "無効なSMILESです。", None, None, None
        
        found_groups = [PKa_RULEBOOK[key] for key, pattern in COMPILED_PATTERNS.items() if pattern and mol.HasSubstructMatch(pattern)]
        
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 300); drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol); drawer.FinishDrawing()
        png_data = drawer.GetDrawingText(); current_png_data = png_data
        img_base_64 = base64.b64encode(png_data).decode('utf-8')

        if not found_groups: return "主要な酸性/塩基性官能基は見つかりませんでした。", img_base_64, "中性に近いと予測されます。", smiles
        
        strongest_acid = min([g for g in found_groups if g["type"] == "acid"], key=lambda x: x["pka"], default=None)
        strongest_base = max([g for g in found_groups if g["type"] == "base"], key=lambda x: x["pka"], default=None)
        
        result_text, details_text = "予測される主要な性質:\n", ""
        if strongest_acid:
            result_text += f"酸性サイト: {strongest_acid['name']} (予測pKa ≈ {strongest_acid['pka']})\n"
            current_predicted_pka = strongest_acid['pka']
        if strongest_base:
            details_text += f"({strongest_base['name']}は塩基性はほぼ示さない)" if strongest_base['name'] == 'アミド' else f"塩基性サイト: {strongest_base['name']} (共役酸の予測pKa ≈ {strongest_base['pka']})"
            if not strongest_acid: current_predicted_pka = strongest_base['pka']
        return result_text.strip(), img_base_64, details_text, smiles

    def calculate_error_click(e):
        try:
            exp_pka = float(exp_pka_input.value)
            if current_predicted_pka is not None:
                error = abs(exp_pka - current_predicted_pka)
                error_display.value = f"文献値との誤差: {error:.2f}"
            else: error_display.value = "先にpKaを予測してください。"
        except: error_display.value = "数値を入力してください。"
        page.update()
    
    def add_to_history(compound):
        user_data = load_user_data()
        if "history" not in user_data: user_data["history"] = []
        history_entry = { "name": compound.iupac_name or (compound.synonyms[0] if compound.synonyms else "不明"), "smiles": compound.smiles }
        user_data["history"] = [item for item in user_data["history"] if item["smiles"] != history_entry["smiles"]]
        user_data["history"].insert(0, history_entry)
        user_data["history"] = user_data["history"][:20]
        save_user_data(user_data)
        update_history_view()

    def update_history_view():
        history_list_view.controls.clear()
        history = load_user_data().get("history", [])
        if not history: history_list_view.controls.append(ft.Text("履歴はありません。"))
        else:
            for item in history: history_list_view.controls.append(ft.ListTile(title=ft.Text(item["name"], size=14), subtitle=ft.Text(item["smiles"], size=10, overflow="ellipsis"), on_click=history_item_click, data=item))
        if page.client_storage: page.update()

    def history_item_click(e):
        smiles = e.control.data["smiles"]
        id_type_dropdown.value = "smiles"; id_textbox.value = smiles
        id_textbox.visible, predict_button.visible, sketcher_view.visible = True, True, False
        run_prediction(smiles)

    def run_prediction(smiles: str):
        compounds = pcp.get_compounds(smiles, 'smiles')
        if not compounds:
            result_display.value, details_display.value, info_container.visible, mol_image.visible = "SMILESから化合物を特定できませんでした。", "", False, False
            page.update(); return
        
        compound = compounds[0]
        add_to_history(compound)
        
        # 新しいプロパティ取得関数を呼び出す
        properties = get_physicochemical_data(compound)

        result_text, img_data, details, retrieved_smiles = predict_pka(smiles)
        result_display.value, details_display.value = result_text, details
        iupac_name_text.value = f"IUPAC名: {compound.iupac_name or 'N/A'}"
        formula_text.value = f"分子式: {compound.molecular_formula or 'N/A'}"
        weight_text.value = f"分子量: {compound.molecular_weight or 'N/A'}"
        cas_text.value = f"CAS番号: {find_cas_number(compound.synonyms)}"
        boiling_point_text.value = f"沸点: {properties.get('boiling_point', 'データなし')}"
        melting_point_text.value = f"融点: {properties.get('melting_point', 'データなし')}"
        density_text.value = f"密度: {properties.get('density', 'データなし')}"
        smiles_text.value = f"SMILES: {retrieved_smiles}"; info_container.visible = True
        
        if img_data: mol_image.src_base64, mol_image.visible = img_data, True
        else: mol_image.visible = False
        
        exp_pka_input.value, error_display.value, new_bp_input.value = "", "", ""
        page.update()

    def add_data_click(e):
        smiles = smiles_text.value.replace("SMILES: ", ""); bp = new_bp_input.value
        if smiles and bp:
            user_data = load_user_data()
            if smiles not in user_data: user_data[smiles] = {}
            if "properties" not in user_data[smiles]: user_data[smiles]["properties"] = {}
            user_data[smiles]["properties"]["boiling_point"] = f"{bp} °C (ユーザー)"
            save_user_data(user_data)
            boiling_point_text.value = f"沸点: {user_data[smiles]['properties']['boiling_point']}"
            page.snack_bar = ft.SnackBar(ft.Text("沸点データを追加しました。"), duration=3000)
            page.snack_bar.open = True; page.update()

    def handle_predict_click(e):
        identifier = id_textbox.value.strip()
        if not identifier:
            result_display.value, details_display.value, info_container.visible, mol_image.visible = "IDを入力してください。", "", False, False
            page.update(); return
        if identifier.upper() == "DCM": run_prediction("C(Cl)Cl"); return
        
        compounds, error = resolve_identifier_to_smiles(identifier, id_type_dropdown.value)
        if error or not compounds:
            result_display.value, details_display.value, info_container.visible, mol_image.visible = error or f"'{identifier}' は見つかりませんでした。", "", False, False
            page.update()
        elif len(compounds) == 1:
            run_prediction(compounds[0].smiles)
        else:
            candidate_dialog.content.controls.clear()
            for comp in compounds[:10]:
                candidate_dialog.content.controls.append(ft.ListTile(title=ft.Text(comp.iupac_name or "名称不明"), subtitle=ft.Text(comp.molecular_formula or "", size=10), on_click=on_candidate_click, data=comp.smiles))
            page.dialog, candidate_dialog.open = candidate_dialog, True; page.update()

    def on_web_message(e): run_prediction(e.data)

    def on_dropdown_change(e):
        result_display.value, details_display.value, info_container.visible, mol_image.visible, id_textbox.value = "", "", False, False, ""
        exp_pka_input.value, error_display.value, new_bp_input.value = "", "", ""
        is_sketcher_mode = (id_type_dropdown.value == "sketch")
        id_textbox.visible, predict_button.visible, sketcher_view.visible = not is_sketcher_mode, not is_sketcher_mode, is_sketcher_mode
        page.update()

    # --- UI Components ---
    title_text = ft.Text("pKa予測アプリ", size=24, weight=ft.FontWeight.BOLD)
    id_type_dropdown = ft.Dropdown(label="入力タイプ", width=400, options=[ft.dropdown.Option("name", "化合物名・略称"), ft.dropdown.Option("smiles", "SMILES"), ft.dropdown.Option("cas", "CAS番号"), ft.dropdown.Option("sketch", "構造を描画")], value="name", on_change=on_dropdown_change)
    id_textbox = ft.TextField(label="IDを入力", hint_text="例: DCM, Aspirin, 50-78-2", width=400, autofocus=True, on_submit=handle_predict_click)
    sketcher_view = WebView("sketcher.html", visible=False, width=420, height=450); sketcher_view.on_web_message = on_web_message
    predict_button = ft.ElevatedButton(text="予測", on_click=handle_predict_click, icon="science", width=200)
    result_display, details_display = ft.Text(size=16, weight=ft.FontWeight.W_500), ft.Text(size=14, color="grey_600")
    
    iupac_name_text, formula_text, weight_text, cas_text, smiles_text = ft.Text(selectable=True), ft.Text(selectable=True), ft.Text(selectable=True), ft.Text(selectable=True), ft.Text(selectable=True)
    boiling_point_text, melting_point_text, density_text = ft.Text(selectable=True), ft.Text(selectable=True), ft.Text(selectable=True)
    
    new_bp_input = ft.TextField(label="沸点を手入力 (°C)", width=200); add_bp_button = ft.ElevatedButton("データを追加", on_click=add_data_click)
    add_data_row = ft.Row([new_bp_input, add_bp_button], alignment=ft.MainAxisAlignment.CENTER)
    exp_pka_input = ft.TextField(label="pKa文献値を入力", width=150); check_error_button = ft.ElevatedButton("誤差を計算", on_click=calculate_error_click)
    error_display = ft.Text(weight=ft.FontWeight.BOLD, selectable=True); error_check_row = ft.Row([exp_pka_input, check_error_button], alignment=ft.MainAxisAlignment.CENTER)
    
    info_container = ft.Container(
        content=ft.Column([
            ft.Text("化合物情報:", weight=ft.FontWeight.BOLD),
            iupac_name_text, formula_text, weight_text, cas_text, 
            boiling_point_text, melting_point_text, density_text, # UIに追加
            smiles_text,
            ft.Divider(), add_data_row, ft.Divider(), error_check_row, error_display,
        ]),
        padding=10, border=ft.border.all(1, "grey_300"), border_radius=5, visible=False, width=400
    )
    mol_image = ft.Image(visible=False, width=300, height=300, fit=ft.ImageFit.CONTAIN); image_container = ft.GestureDetector(content=mol_image, on_secondary_tap_down=open_context_menu)
    history_list_view = ft.ListView(expand=True, spacing=5, padding=ft.padding.only(top=10))
    history_container = ft.Container(content=ft.Column([ft.Text("検索履歴", weight=ft.FontWeight.BOLD), history_list_view]), expand=True, padding=10, border=ft.border.all(1, "outline"), border_radius=5)

    main_column = ft.Column(
        [
            ft.Container(height=10), title_text, ft.Container(height=10),
            id_type_dropdown, id_textbox, sketcher_view, predict_button,
            ft.Container(height=10), ft.Divider(),
            result_display, details_display, ft.Container(height=10),
            info_container, image_container,
        ],
        scroll=ft.ScrollMode.ADAPTIVE, expand=True, horizontal_alignment=ft.CrossAxisAlignment.CENTER, spacing=10,
    )
    page.add(ft.Row([ft.Container(content=history_container, width=280, padding=10), ft.VerticalDivider(), ft.Container(content=main_column, expand=True, padding=10)], expand=True))
    update_history_view()

if __name__ == "__main__":
    ft.app(target=main)