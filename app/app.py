"""
Cruzipain-IC50Pred - Aplicacao de triagem virtual de inibidores de Cruzipain (T. cruzi)

Pipeline de predicao em tres etapas:
  1. Triagem    -> classifica a molecula como inibidor ou decoy (Modelo 1)
  2. Potencia   -> classifica o inibidor em Classe 0/1/2 (Modelo 2)
  3. Regressao  -> estima o valor de IC50 em nM (Modelo 3, via pIC50 + rank)

Atalho de molecula conhecida:
  Se o SMILES informado ja existe no conjunto de treino, o resultado e
  recuperado diretamente da base, sem passar pelos modelos.

Como executar:
  streamlit run app/app.py
"""

from pathlib import Path
import io

import numpy as np
import pandas as pd
import pickle
import streamlit as st
import plotly.express as px
from padelpy import from_smiles


# --------------------------------------------------------------------------
# Configuracao de caminhos
# --------------------------------------------------------------------------
BASE = Path(__file__).resolve().parent.parent
MODELS = BASE / "models"
DATA = BASE / "data" / "processed"


# --------------------------------------------------------------------------
# Carregamento de modelos e referencias (cacheado)
# --------------------------------------------------------------------------
@st.cache_resource
def load_artifacts():
    """Carrega os tres modelos, transformadores e a base de moleculas conhecidas."""
    art = {}

    with open(MODELS / "decoy_inhibitor_rf.pkl", "rb") as f:
        art["model1"] = pickle.load(f)
    with open(MODELS / "imputer_model.pkl", "rb") as f:
        art["imputer"] = pickle.load(f)
    with open(MODELS / "selected_features.pkl", "rb") as f:
        art["features"] = pickle.load(f)            # 628 descritores
    with open(MODELS / "HGB_model_potency.pkl", "rb") as f:
        art["model2"] = pickle.load(f)
    with open(MODELS / "rf_model_regression.pkl", "rb") as f:
        art["model3"] = pickle.load(f)
    with open(MODELS / "regression_coefficients.pkl", "rb") as f:
        art["reg_coef"] = pickle.load(f)            # polinomio rank -> pIC50
    with open(MODELS / "regression_features.pkl", "rb") as f:
        art["reg_features"] = pickle.load(f)        # 628 features + 'Class'

    # Base de moleculas conhecidas (treino) - usada no atalho e na referencia de rank
    known = pd.read_csv(DATA / "inhibitors_pruned.csv")
    art["known"] = known

    # Matriz de referencia para o ranqueamento da regressao (628 features)
    art["rank_ref"] = known[art["features"]].apply(pd.to_numeric, errors="coerce")
    art["rank_ref"] = (art["rank_ref"]
                       .replace([np.inf, -np.inf], np.nan)
                       .clip(-1e6, 1e6))
    art["rank_ref"] = art["rank_ref"].fillna(art["rank_ref"].median())

    # Medianas de treino para imputacao de moleculas novas
    art["medians"] = art["rank_ref"].median()

    # Limites de capping IQR fixos, salvos no treino do Modelo 1.
    # Treino e producao usam exatamente os mesmos limites (sem train/serving skew).
    with open(MODELS / "cap_limits.pkl", "rb") as f:
        cap = pickle.load(f)
    art["cap_lower"] = cap["lower"]
    art["cap_upper"] = cap["upper"]
    if "medians" in cap:
        art["cap_medians"] = cap["medians"]
    else:
        art["cap_medians"] = art["medians"]

    return art


def normalize_smiles(smiles: str) -> str:
    """Normaliza o SMILES para comparacao textual (remove espacos).
    Retorna string vazia se for vazio/invalido, para nao gerar match falso."""
    s = (smiles or "").strip()
    if not s or s.lower() == "nan":
        return ""
    return s


# --------------------------------------------------------------------------
# Etapas do pipeline
# --------------------------------------------------------------------------
def smiles_to_descriptors(smiles: str) -> pd.Series:
    """Calcula os 1875 descritores PaDEL (com 3D) para um SMILES."""
    desc = from_smiles(smiles, timeout=120)
    return pd.Series(desc, dtype="object")


def align_and_clean(desc: pd.Series, art) -> pd.DataFrame:
    """Alinha aos 628 descritores selecionados e trata valores numericos."""
    row = pd.DataFrame([desc]).apply(pd.to_numeric, errors="coerce")
    row = row.reindex(columns=art["features"])
    row = row.replace([np.inf, -np.inf], np.nan)
    # Preenche faltantes com a mediana de treino
    row = row.fillna(art["medians"])
    row = row.clip(-1e6, 1e6)
    return row


def prep_model1(row: pd.DataFrame, art) -> pd.DataFrame:
    """Pre-processamento identico ao treino do Modelo 1: capping fixo + median + imputer."""
    capped = row.clip(lower=art["cap_lower"], upper=art["cap_upper"], axis=1)
    capped = capped.replace([np.inf, -np.inf], np.nan).fillna(art["cap_medians"]).clip(-1e6, 1e6)
    return pd.DataFrame(art["imputer"].transform(capped), columns=capped.columns)


def predict_type(row: pd.DataFrame, art) -> int:
    """Modelo 1: 1 = inibidor, 0 = decoy.
    Usa o mesmo pre-processamento do treino (capping fixo + imputer)."""
    x = prep_model1(row, art)
    return int(art["model1"].predict(x)[0])


def predict_class(row: pd.DataFrame, art) -> int:
    """Modelo 2: classe de potencia 0, 1 ou 2."""
    return int(art["model2"].predict(row)[0])


def predict_ic50(row: pd.DataFrame, potency_class: int, art) -> float:
    """Modelo 3: estima IC50 (nM) via rank do pIC50 e polinomio de conversao."""
    ref = art["rank_ref"]
    n = len(ref)

    # Ranqueia cada feature da molecula nova em relacao a base de treino
    ranked = {}
    for feat in art["features"]:
        v = row.iloc[0][feat]
        ranked[feat] = int((ref[feat].values < v).sum()) + 1

    # A classe entra como atributo (mesmo papel do treino)
    ranked["Class"] = potency_class

    x = pd.DataFrame([ranked])[art["reg_features"]]

    # rank previsto -> pIC50 -> IC50 (nM)
    rank_pred = art["model3"].predict(x)[0]
    pic50 = np.polyval(art["reg_coef"], rank_pred)
    ic50 = 10 ** (9 - pic50)
    return float(ic50)


def lookup_known(smiles: str, art):
    """Procura a molecula na base de treino. Retorna a linha ou None."""
    target = normalize_smiles(smiles)
    if not target:                      # SMILES vazio nao casa com nada
        return None
    known = art["known"]
    for _, r in known.iterrows():
        ref = normalize_smiles(str(r["smiles"]))
        if ref and ref == target:       # so casa se ambos forem validos e iguais
            return r
    return None


# --------------------------------------------------------------------------
# Processamento de uma lista de moleculas
# --------------------------------------------------------------------------
def run_pipeline(entries, art):
    """entries: lista de (smiles, nome). Retorna DataFrame de resultados."""
    results = []
    prog = st.progress(0.0, text="Processando moleculas...")

    for i, (smi, name) in enumerate(entries):
        item = {"ID": i + 1, "Composto": smi, "Nome": name or "N/A"}

        # 1) Atalho: molecula conhecida
        hit = lookup_known(smi, art)
        if hit is not None:
            cls = int(hit["potency_class"])
            item.update({
                "Tipo": "Inibidor",
                "Classe": cls,
                "IC50 (nM)": round(float(hit["IC50_nM"]), 2),
                "Origem": "Base conhecida",
            })
            results.append(item)
            prog.progress((i + 1) / len(entries))
            continue

        # 2) Molecula nova: pipeline completo
        try:
            desc = smiles_to_descriptors(smi)
            row = align_and_clean(desc, art)
        except Exception as e:
            item.update({"Tipo": "Erro", "Classe": "-",
                         "IC50 (nM)": "-", "Origem": f"{type(e).__name__}: {e}"})
            results.append(item)
            prog.progress((i + 1) / len(entries))
            continue

        tipo = predict_type(row, art)
        if tipo == 0:
            item.update({"Tipo": "Decoy", "Classe": "N/A",
                         "IC50 (nM)": "N/A", "Origem": "Modelo"})
        else:
            cls = predict_class(row, art)
            ic50 = predict_ic50(row, cls, art)
            item.update({"Tipo": "Inibidor", "Classe": cls,
                         "IC50 (nM)": round(ic50, 2), "Origem": "Modelo"})

        results.append(item)
        prog.progress((i + 1) / len(entries))

    prog.empty()
    return pd.DataFrame(results)


# --------------------------------------------------------------------------
# Parsing da entrada
# --------------------------------------------------------------------------
def parse_text(text: str):
    """Le 'SMILES, Nome' por linha. Maximo de 20 compostos."""
    entries = []
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else ""
        entries.append((smi, name))
    return entries[:20]


def parse_file(uploaded):
    """Le CSV/XLSX: SMILES na 1a coluna, nome na 2a (opcional)."""
    if uploaded.name.endswith((".xlsx", ".xls")):
        df = pd.read_excel(uploaded)
    else:
        df = pd.read_csv(uploaded)
    entries = []
    for _, r in df.iterrows():
        smi = str(r.iloc[0]).strip()
        name = str(r.iloc[1]).strip() if df.shape[1] > 1 else ""
        if smi and smi.lower() != "nan":
            entries.append((smi, name))
    return entries[:20]


# --------------------------------------------------------------------------
# Apresentacao dos resultados
# --------------------------------------------------------------------------
CLASS_LABELS = {
    0: "Classe 0 - Mais potente",
    1: "Classe 1 - Moderado",
    2: "Classe 2 - Fraco",
}
CLASS_COLORS = {0: "#1D9E75", 1: "#8E6FD8", 2: "#E24B82"}


def show_results(df):
    st.markdown("## Resumo dos resultados")

    # Legenda das classes de potencia
    st.markdown(
        "**Classes de potencia dos inibidores:** "
        ":green[Classe 0 - Mais potente] &nbsp;|&nbsp; "
        ":violet[Classe 1 - Moderado] &nbsp;|&nbsp; "
        ":red[Classe 2 - Fraco]"
    )

    col1, col2 = st.columns(2)

    # Distribuicao de tipo (inibidor vs decoy)
    with col1:
        tipo_counts = df["Tipo"].value_counts()
        fig_pie = px.pie(
            names=tipo_counts.index, values=tipo_counts.values,
            title="Distribuicao por tipo",
            color=tipo_counts.index,
            color_discrete_map={"Inibidor": "#F5A623", "Decoy": "#4A90D9",
                                "Erro": "#999999"},
            hole=0.45,
        )
        fig_pie.update_layout(height=330, margin=dict(t=50, b=10, l=10, r=10))
        st.plotly_chart(fig_pie, use_container_width=True)

    # Distribuicao das classes de inibidores
    with col2:
        inh = df[df["Tipo"] == "Inibidor"].copy()
        if not inh.empty:
            inh = inh[pd.to_numeric(inh["Classe"], errors="coerce").notna()]
            cls_counts = inh["Classe"].astype(int).value_counts().sort_index()
            fig_bar = px.bar(
                x=[CLASS_LABELS[c] for c in cls_counts.index],
                y=cls_counts.values,
                title="Distribuicao das classes de potencia",
                color=[CLASS_LABELS[c] for c in cls_counts.index],
                color_discrete_map={CLASS_LABELS[c]: CLASS_COLORS[c]
                                    for c in cls_counts.index},
            )
            fig_bar.update_layout(height=330, showlegend=False,
                                  margin=dict(t=50, b=10, l=10, r=10),
                                  xaxis_title="", yaxis_title="Contagem")
            st.plotly_chart(fig_bar, use_container_width=True)
        else:
            st.info("Nenhum inibidor identificado nesta submissao.")

    # Tabela de resultados
    st.markdown("### Detalhamento")
    st.dataframe(df, use_container_width=True, hide_index=True)

    # Exportar CSV
    buffer = io.StringIO()
    df.to_csv(buffer, index=False)
    st.download_button(
        "Exportar resultados em CSV",
        data=buffer.getvalue(),
        file_name="cruzipain_ic50pred_resultados.csv",
        mime="text/csv",
    )


# --------------------------------------------------------------------------
# Interface
# --------------------------------------------------------------------------
def main():
    st.set_page_config(page_title="Cruzipain-IC50Pred", page_icon="🧬",
                       layout="wide")

    st.markdown(
        "<h1 style='text-align:center;color:#7B2D8E;'>Cruzipain-IC50Pred</h1>",
        unsafe_allow_html=True,
    )
    st.markdown(
        "<p style='text-align:center;color:#555;'>Classifica moleculas em "
        "inibidores e decoys, categoriza os inibidores por potencia e estima "
        "o valor de IC50 contra a Cruzipain de <i>Trypanosoma cruzi</i>.</p>",
        unsafe_allow_html=True,
    )

    try:
        art = load_artifacts()
    except Exception as e:
        st.error(f"Falha ao carregar os modelos: {e}")
        st.stop()

    tab_smiles, tab_file = st.tabs(["Inserir SMILES", "Enviar arquivo"])

    entries = None
    with tab_smiles:
        text = st.text_area(
            "SMILES (um por linha, formato: SMILES, Nome)",
            height=180,
            placeholder="CCC, Composto A\n"
                        "CNC(=O)C1=CN=CN1, Composto B\n"
                        "Um composto por linha. Maximo de 20.",
        )
        if st.button("Prever", type="primary", key="btn_smiles"):
            entries = parse_text(text)

    with tab_file:
        uploaded = st.file_uploader(
            "CSV ou XLSX - SMILES na primeira coluna, nome na segunda (opcional)",
            type=["csv", "xlsx", "xls"],
        )
        if st.button("Prever", type="primary", key="btn_file"):
            if uploaded is not None:
                entries = parse_file(uploaded)
            else:
                st.warning("Envie um arquivo antes de prever.")

    if entries is not None:
        if not entries:
            st.warning("Nenhum SMILES valido informado.")
        else:
            df = run_pipeline(entries, art)
            show_results(df)


if __name__ == "__main__":
    main()