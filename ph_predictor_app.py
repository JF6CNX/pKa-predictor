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
        json.dump(data, f, indent=4)

# --- NISTから沸点をスクレイピングする関数 ---
def scrape_nist_boiling_point(cas_number):
    if not cas_number or cas_number == "N/A": return None
    try:
        url = f"https://webbook.nist.gov/cgi/cbook.cgi?ID={cas_number}&Units=SI"
        response = requests.get(url, timeout=5)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        bp_header = soup.find('th', string=lambda t: t and 'Boiling Point' in t)
        if bp_header:
            bp_value_td = bp_header.find_next_sibling('td')
            if bp_value_td: return f"{bp_value_td.get_text(strip=True)} °C (NIST)"
        return None
    except Exception:
        return None


def main(page: ft.Page):
    page.title = "pKa Predictor with History"
    page.window_width = 800
    page.window_height = 900
    page.theme_mode = ft.ThemeMode.LIGHT
    
    current_png_data = None
    current_predicted_pka = None

    def save_file_result(e: ft.FilePickerResultEvent):
        if e.path and current_png_data:
            with open(e.path, "wb") as f: f.write(current_png_data)
            page.snack_bar = ft.SnackBar(ft.Text(f"画像を {e.path} に保存しました。"), duration=3000)
            page.snack_bar.open = True; page.update()
    save_dialog = ft.FilePicker(on_result=save_file_result)
    
    # --- 右クリックメニュー関連の関数 ---
    def open_context_menu(e: ft.TapEvent):
        if current_png_data and mol_image.visible:
            context_menu.left = e.global_x
            context_menu.top = e.global_y
            context_menu.visible = True
            page.update()

    def close_context_menu(e):
        if context_menu.visible:
            context_menu.visible = False
            page.update()

    def save_image_click(e):
        close_context_menu(e)
        if current_png_data:
            save_dialog.save_file(
                dialog_title="画像を保存",
                file_name=f"{(id_textbox.value or 'molecule').replace(' ', '_')}.png",
                allowed_extensions=["png"]
            )
        else:
            page.snack_bar = ft.SnackBar(ft.Text("保存する画像がありません。先に構造を予測してください。"), duration=3000, bgcolor="orange_accent_700")
            page.snack_bar.open = True
            page.update()
    
    context_menu = ft.Card(
        content=ft.Column([
            ft.ListTile(
                title=ft.Text("画像を保存する", color="black"),
                on_click=save_image_click,
                autofocus=True,
                hover_color="grey_200",
                width=200,
            )
        ], spacing=0),
        elevation=4,
        visible=False,
        color="white",
        shape=ft.RoundedRectangleBorder(radius=ft.border_radius.all(8)),
    )
    
    page.overlay.extend([save_dialog, context_menu])
    page.on_click = close_context_menu

    def find_cas_number(synonyms):
        cas_pattern = re.compile(r'^\d{2,7}-\d{2}-\d$')
        for syn in synonyms:
            if cas_pattern.match(syn): return syn
        return "N/A"
    def get_compound_from_identifier(identifier, id_type):
        try:
            results = pcp.get_compounds(identifier, id_type)
            if not results: return None, f"{id_type} '{identifier}' は見つかりませんでした。"
            return results[0], None
        except Exception as e:
            return None, f"検索中にエラーが発生しました: {e}"

    def get_boiling_point(compound):
        user_data = load_user_data()
        smiles = compound.connectivity_smiles
        if smiles in user_data and "boiling_point" in user_data[smiles]:
            return user_data[smiles]["boiling_point"]
        try:
            bp = compound.boiling_point
            if bp: return f"{bp} °C (PubChem)"
        except (KeyError, AttributeError): pass
        try:
            for prop in compound.to_dict(properties=['synonyms']).get('synonyms', []):
                if 'Boiling Point' in prop:
                    value = prop.split(':')[-1].strip()
                    if value: return value
        except Exception: pass
        cas = find_cas_number(compound.synonyms)
        nist_bp = scrape_nist_boiling_point(cas)
        if nist_bp: return nist_bp
        return "N/A"

    def predict_pka(smiles: str):
        nonlocal current_png_data, current_predicted_pka
        current_predicted_pka = None
        current_png_data = None
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: raise ValueError("Invalid SMILES")
        except: return "無効なSMILESです。", None, None, None
        
        found_groups = []
        for key, pattern in COMPILED_PATTERNS.items():
            if pattern and mol.HasSubstructMatch(pattern):
                found_groups.append(PKa_RULEBOOK[key])
        
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
        drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol); drawer.FinishDrawing()
        png_data = drawer.GetDrawingText(); current_png_data = png_data
        
        temp_png_path = "temp_mol.png"
        with open(temp_png_path, "wb") as f: f.write(png_data)
        with open(temp_png_path, "rb") as f: img_base_64 = base64.b64encode(f.read()).decode('utf-8')
        os.remove(temp_png_path)

        if not found_groups: return "主要な酸性/塩基性官能基は見つかりませんでした。", img_base_64, "中性に近いと予測されます。", smiles
        strongest_acid = min([g for g in found_groups if g["type"] == "acid"], key=lambda x: x["pka"], default=None)
        strongest_base = max([g for g in found_groups if g["type"] == "base"], key=lambda x: x["pka"], default=None)
        result_text = "予測される主要な性質:\n"; details_text = ""
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
            else:
                error_display.value = "先にpKaを予測してください。"
        except (ValueError, TypeError):
            error_display.value = "数値を入力してください。"
        page.update()
        
    def add_to_history(compound):
        user_data = load_user_data()
        if "history" not in user_data: user_data["history"] = []
        history_entry = { "name": compound.iupac_name or (compound.synonyms[0] if compound.synonyms else "Unknown"), "smiles": compound.connectivity_smiles }
        user_data["history"] = [item for item in user_data["history"] if item["smiles"] != history_entry["smiles"]]
        user_data["history"].insert(0, history_entry)
        user_data["history"] = user_data["history"][:20]
        save_user_data(user_data)
        update_history_view()

    def update_history_view():
        history_list_view.controls.clear()
        user_data = load_user_data()
        history = user_data.get("history", [])
        if not history:
            history_list_view.controls.append(ft.Text("履歴はありません。"))
        else:
            for item in history:
                history_list_view.controls.append(ft.ListTile(title=ft.Text(item["name"], size=14),
                    subtitle=ft.Text(item["smiles"], size=10, overflow="ellipsis"), on_click=history_item_click, data=item))
        if page.client_storage: page.update()

    def history_item_click(e):
        history_item = e.control.data
        smiles = history_item["smiles"]; name = history_item["name"]
        id_type_dropdown.value = "name"; id_textbox.value = name
        id_textbox.visible = True; predict_button.visible = True; sketcher_view.visible = False
        run_prediction(smiles)

    def run_prediction(smiles: str):
        compound, error = get_compound_from_identifier(smiles, "smiles")
        if error:
            result_display.value = error; details_display.value = ""; info_container.visible = False; mol_image.visible = False
            page.update(); return
        
        add_to_history(compound)
        
        result_text, img_data, details, retrieved_smiles = predict_pka(smiles)
        result_display.value = result_text; details_display.value = details
        iupac_name_text.value = f"IUPAC名: {compound.iupac_name or 'N/A'}"
        formula_text.value = f"分子式: {compound.molecular_formula or 'N/A'}"
        weight_text.value = f"分子量: {compound.molecular_weight or 'N/A'}"
        cas_text.value = f"CAS番号: {find_cas_number(compound.synonyms)}"
        boiling_point_text.value = f"沸点: {get_boiling_point(compound)}"
        smiles_text.value = f"SMILES: {retrieved_smiles}"; info_container.visible = True
        if img_data:
            mol_image.src_base64 = img_data; mol_image.visible = True
        else: mol_image.visible = False
        exp_pka_input.value = ""; error_display.value = ""
        new_bp_input.value = ""
        page.update()

    def add_data_click(e):
        smiles = smiles_text.value.replace("SMILES: ", "")
        bp = new_bp_input.value
        if smiles and bp:
            user_data = load_user_data()
            if smiles not in user_data: user_data[smiles] = {}
            user_data[smiles]["boiling_point"] = f"{bp} °C (user)"
            save_user_data(user_data)
            boiling_point_text.value = f"沸点: {user_data[smiles]['boiling_point']}"
            page.snack_bar = ft.SnackBar(ft.Text("沸点データを追加しました。"), duration=3000)
            page.snack_bar.open = True; page.update()

    def handle_predict_click(e):
        identifier = id_textbox.value
        id_type = id_type_dropdown.value
        compound, error = get_compound_from_identifier(identifier, id_type)
        if error:
            result_display.value = error; details_display.value = ""; info_container.visible = False; mol_image.visible = False
            page.update()
        else: run_prediction(compound.connectivity_smiles)

    def on_web_message(e):
        smiles = e.data
        run_prediction(smiles)

    def on_dropdown_change(e):
        result_display.value = ""; details_display.value = ""; info_container.visible = False
        mol_image.visible = False; id_textbox.value = ""
        exp_pka_input.value = ""; error_display.value = ""
        new_bp_input.value = ""
        is_sketcher_mode = (id_type_dropdown.value == "sketch")
        id_textbox.visible = not is_sketcher_mode
        predict_button.visible = not is_sketcher_mode
        sketcher_view.visible = is_sketcher_mode
        page.update()

    # --- UIコンポーネント定義 ---
    title_text = ft.Text("pKa Predictor", size=24, weight=ft.FontWeight.BOLD)
    id_type_dropdown = ft.Dropdown(label="入力タイプ", width=400, options=[
        ft.dropdown.Option("name", "化合物名 (英語)"), ft.dropdown.Option("smiles", "SMILES"),
        ft.dropdown.Option("cas", "CAS番号"), ft.dropdown.Option("sketch", "構造を描画"),
    ], value="name", on_change=on_dropdown_change)
    id_textbox = ft.TextField(label="IDを入力", hint_text="例: Aspirin, 50-78-2", width=400, autofocus=True, on_submit=handle_predict_click)
    sketcher_view = WebView("sketcher.html", visible=False, width=420, height=450); sketcher_view.on_web_message = on_web_message
    predict_button = ft.ElevatedButton(text="予測する", on_click=handle_predict_click, icon="science", width=200)
    result_display, details_display = ft.Text(size=16, weight=ft.FontWeight.W_500), ft.Text(size=14, color="grey_600")
    iupac_name_text, formula_text, weight_text, cas_text, smiles_text, boiling_point_text = ft.Text(), ft.Text(), ft.Text(), ft.Text(), ft.Text(selectable=True), ft.Text()
    new_bp_input = ft.TextField(label="沸点を手入力 (°C)", width=200); add_bp_button = ft.ElevatedButton("データを追加", on_click=add_data_click)
    add_data_row = ft.Row([new_bp_input, add_bp_button], alignment=ft.MainAxisAlignment.CENTER)
    exp_pka_input = ft.TextField(label="pKa文献値を入力", width=150)
    check_error_button = ft.ElevatedButton("誤差を表示", on_click=calculate_error_click)
    error_display = ft.Text(weight=ft.FontWeight.BOLD)
    error_check_row = ft.Row([exp_pka_input, check_error_button], alignment=ft.MainAxisAlignment.CENTER)
    info_container = ft.Container(
        content=ft.Column([
            ft.Text("化合物情報:", weight=ft.FontWeight.BOLD),
            iupac_name_text, formula_text, weight_text, cas_text, boiling_point_text, smiles_text,
            ft.Divider(), add_data_row, ft.Divider(), error_check_row, error_display,
        ]),
        padding=10, border=ft.border.all(1, "grey_300"), border_radius=5, visible=False, width=400)
    mol_image = ft.Image(visible=False, width=300, height=300, fit=ft.ImageFit.CONTAIN)
    
    image_container = ft.GestureDetector(
        content=mol_image,
        on_secondary_tap_down=open_context_menu
    )
        
    history_list_view = ft.ListView(expand=True, spacing=5, padding=ft.padding.only(top=10))
    history_container = ft.Container(
        content=ft.Column([ft.Text("検索履歴", weight=ft.FontWeight.BOLD), history_list_view]),
        expand=True, padding=10, border=ft.border.all(1, "outline"), border_radius=5,
    )

    # --- UIレイアウト ---
    main_column = ft.Column(
        [
            ft.Container(height=10), title_text, ft.Container(height=10),
            id_type_dropdown, id_textbox, sketcher_view, predict_button,
            ft.Container(height=10), ft.Divider(),
            result_display, details_display, ft.Container(height=10),
            info_container, image_container,
        ],
        scroll=ft.ScrollMode.ADAPTIVE, expand=True,
        horizontal_alignment=ft.CrossAxisAlignment.CENTER, spacing=10,
    )
    page.add(ft.Row([
        ft.Container(content=history_container, width=280, padding=10),
        ft.VerticalDivider(),
        ft.Container(content=main_column, expand=True, padding=10),
    ], expand=True))
    update_history_view()

ft.app(target=main)