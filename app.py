import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import numpy as np
import os
import pubchempy as pcp

def name_to_smiles(name):
    try:
        compound = pcp.get_compounds(name, 'name')[0]
        return compound.isomeric_smiles
    except:
        return None

# --- Deskriptoren berechnen ---
def calc_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolWt": Descriptors.MolWt(mol),
        "logP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
        "FractionCSP3": Descriptors.FractionCSP3(mol),
        "RingCount": Descriptors.RingCount(mol),
        "MolMR": Descriptors.MolMR(mol),
        "NumAromaticRings": Descriptors.NumAromaticRings(mol),
    }

# --- RT-Anpassung an Bedingungen ---
def adjust_rt(base_rt, temp, ph, säule):
    rt = base_rt * (0.98 ** (temp - 25))
    if ph < 4:
        rt *= 1.05
    elif ph > 8:
        rt *= 0.95
    if säule == "Phenyl":
        rt *= 1.1
    elif säule == "HILIC":
        rt *= 0.8
    return rt

# --- Rs berechnen ---
def berechne_rs(rt1, rt2, w=0.5):
    return (2 * abs(rt1 - rt2)) / (2 * w)

# --- Chromatogramm zeichnen ---
def plot_peaks(names, rts):
    fig, ax = plt.subplots()
    for name, rt in zip(names, rts):
        x = np.linspace(rt - 1.5, rt + 1.5, 300)
        y = np.exp(-((x - rt)**2) / (2 * 0.25**2))
        ax.plot(x, y, label=name)
    ax.set_xlabel("Zeit (min)")
    ax.set_ylabel("Signal")
    ax.set_title("Simuliertes HPLC-Chromatogramm")
    ax.legend()
    st.pyplot(fig)
# --- Streamlit-Grundlayout ---
st.set_page_config(page_title="HPLC-Tool", layout="wide")
st.title("🔬 Intelligentes HPLC-Peak-Vorhersage-Tool")

st.sidebar.header("🔧 Messbedingungen")
temp = st.sidebar.slider("Temperatur [°C]", 20, 60, 25)
ph = st.sidebar.slider("pH-Wert", 1, 14, 7)
säule = st.sidebar.selectbox("Säulenmaterial", ["C18", "Phenyl", "HILIC"])

# --- Interne Wirkstoffdatenbank ---
substanz_db = {
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Diclofenac": "COC1=CC=C(C=C1)C2=CC=CC=C2Cl",
    "Metoprolol": "CC(C)NCC(COC1=CC=CC=C1)O",
    "Lidocain": "CCN(CC)C(=O)C1=CC=CC=C1N(C)C",
    "Naproxen": "CC(C)C1=CC=C(C=C1)C(O)C(=O)O",
    "Amlodipin": "CCOC(=O)C1=C(C(=C(N=C1C)C)C(=O)OCC)N",
    "Omeprazol": "CC1=C(N=CN=C1)S(=O)C2=CC=CC=C2OC3=CC=CC=C3"
    # Du kannst hier bis zu 50+ erweitern...
}

# --- Auswahl: intern + PubChem + Löschen ---

# SessionState initialisieren
if "pubchem_inputs" not in st.session_state:
    st.session_state.pubchem_inputs = []
    st.session_state.pubchem_names = []

# Auswahl aus interner Datenbank
st.markdown("### 📋 Wähle bekannte Substanzen aus interner Liste:")
auswahl = st.multiselect("Substanzen auswählen:", list(substanz_db.keys()))
interne_inputs = [substanz_db[name] for name in auswahl]
interne_namen = auswahl

# PubChem-Suche
st.markdown("### 🔎 Weitere Substanzen per Namen aus PubChem hinzufügen:")
eingabe_name = st.text_input("Wirkstoffname (einzeln eingeben)")

col1, col2 = st.columns([1, 1])

with col1:
    if st.button("➕ Hinzufügen aus PubChem"):
        gefunden = name_to_smiles(eingabe_name)
        if gefunden:
            st.session_state.pubchem_inputs.append(gefunden)
            st.session_state.pubchem_names.append(eingabe_name)
            st.success(f"Hinzugefügt: {eingabe_name}")
        else:
            st.error("Nicht gefunden.")

with col2:
    if st.button("🗑️ PubChem-Eingaben löschen"):
        st.session_state.pubchem_inputs = []
        st.session_state.pubchem_names = []
        st.info("Alle PubChem-Eingaben gelöscht.")

# Zusammenführen aller Substanzen
valid_inputs = interne_inputs + st.session_state.pubchem_inputs
namen = interne_namen + st.session_state.pubchem_names

# Anzeige der aktuellen Liste
if namen:
    st.markdown("### 📦 Aktuelle Substanzen im Vergleich:")
    for i, (name, smi) in enumerate(zip(namen, valid_inputs), 1):
        st.write(f"{i}. **{name}** – {smi}")

st.markdown("### 🔎 SMILES automatisch per Wirkstoffname laden")
eingabe_name = st.text_input("Wirkstoffname eingeben (z. B. Ibuprofen)")

if st.button("🔄 SMILES von PubChem holen"):
    gefunden = name_to_smiles(eingabe_name)
    if gefunden:
        st.success(f"SMILES für {eingabe_name}: {gefunden}")
        valid_inputs = [gefunden]
        namen = [eingabe_name]
    else:
        st.error("Wirkstoff nicht gefunden.")

# --- Hochladen und Anzeigen einer externen Excel-Datei ---
st.markdown("### 📁 Optional: Retentionszeit-Datenbank (Excel) hochladen")

hochgeladene_datei = st.file_uploader("Lade eine Excel-Datei hoch (mit Spalten: Substanz, SMILES, RT [min], Säule, pH, Temperatur):", type=["xlsx"])

if hochgeladene_datei is not None:
    with open("daten/rt_daten.xlsx", "wb") as f:
        f.write(hochgeladene_datei.read())
    st.success("Datei gespeichert! Starte die App neu, um neue Daten zu laden.")

# --- Laden vorhandener Excel-Datei (wenn vorhanden) ---
dateipfad = "daten/rt_daten.xlsx"
if os.path.exists(dateipfad):
    try:
        df_rt = pd.read_excel(dateipfad)
        st.markdown("### 📊 Geladene RT-Daten (Excel-Datei)")
        st.dataframe(df_rt)
    except:
        st.warning("Fehler beim Laden der Datei. Ist sie korrekt formatiert?")
else:
    st.info("Noch keine RT-Datei gefunden.")
# --- Vorhersage starten ---
if st.button("🔍 Retentionszeiten & Trennung berechnen"):
    if not valid_inputs:
        st.warning("Bitte mindestens eine Substanz auswählen.")
    else:
        results = []
        for smi in valid_inputs:
            desc = calc_descriptors(smi)
            if desc is None:
                st.error(f"Ungültige SMILES: {smi}")
                continue
            base_rt = 0.03 * desc["MolWt"] + 0.8 * desc["logP"]
            rt = adjust_rt(base_rt, temp, ph, säule)
            results.append({"SMILES": smi, "RT": rt})

        df = pd.DataFrame(results)
        df["Substanz"] = namen

        st.subheader("📋 Vorhergesagte Retentionszeiten")
        st.dataframe(df[["Substanz", "SMILES", "RT"]].round(2))

        st.subheader("🧪 Trennbarkeitsbewertung")
        for i in range(len(df)):
            for j in range(i+1, len(df)):
                rs = berechne_rs(df.loc[i, "RT"], df.loc[j, "RT"])
                status = "✅ gut getrennt" if rs >= 1.5 else ("⚠️ knapp" if rs >= 1.0 else "❌ überlappt")
                st.write(f"{df.loc[i, 'Substanz']} vs. {df.loc[j, 'Substanz']} → Rs = {rs:.2f} {status}")

        st.subheader("📈 Simuliertes Chromatogramm")
        plot_peaks(df["Substanz"], df["RT"])
# --- Optimale Bedingungen finden ---
def finde_optimale_bedingungen(smiles_liste):
    best_combo = None
    best_min_rs = -1
    best_rts = []
    temperaturen = range(20, 65, 5)
    ph_werte = range(3, 11)
    säulen = ["C18", "Phenyl", "HILIC"]

    for t in temperaturen:
        for ph in ph_werte:
            for säule in säulen:
                rts = []
                for smi in smiles_liste:
                    desc = calc_descriptors(smi)
                    if desc is None:
                        continue
                    base_rt = 0.03 * desc["MolWt"] + 0.8 * desc["logP"]
                    rt = adjust_rt(base_rt, t, ph, säule)
                    rts.append(rt)
                if len(rts) < 2:
                    continue
                min_rs = min([
                    berechne_rs(rts[i], rts[j])
                    for i in range(len(rts)) for j in range(i+1, len(rts))
                ])
                if min_rs > best_min_rs:
                    best_min_rs = min_rs
                    best_combo = (t, ph, säule)
                    best_rts = rts
    return best_combo, best_min_rs, best_rts

# --- Button zur automatischen Optimierung ---
st.markdown("---")
st.markdown("### 🧠 Automatische Optimierung der Messbedingungen")
if st.button("⚙️ Beste Bedingungen berechnen"):
    if not valid_inputs:
        st.warning("Bitte zuerst Substanzen auswählen.")
    else:
        optimal, min_rs, rts_opt = finde_optimale_bedingungen(valid_inputs)
        if optimal:
            st.success(f"Optimale Trennung bei: Temperatur = {optimal[0]} °C, pH = {optimal[1]}, Säule = {optimal[2]}")
            st.write(f"Minimaler Rs unter diesen Bedingungen: {min_rs:.2f}")
            plot_peaks(namen, rts_opt)
        else:
            st.error("Optimierung fehlgeschlagen – nicht genug gültige Substanzen.")

