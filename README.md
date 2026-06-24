# Chagas Drug Discovery — Triagem Virtual de Inibidores de Cruzipaína

Pipeline de aprendizado de máquina para a descoberta de inibidores da **Cruzipaína** (Cruzain), a principal cisteíno-protease do *Trypanosoma cruzi*, agente causador da doença de Chagas. O projeto reproduz e adapta uma metodologia de triagem virtual em três etapas, originalmente desenvolvida para um alvo de Alzheimer, aplicando-a a um novo alvo terapêutico.

A doença de Chagas afeta milhões de pessoas, sobretudo na América Latina, e dispõe de poucas opções terapêuticas. A Cruzipaína é essencial para a sobrevivência e a replicação do parasita, o que a torna um alvo molecular promissor.

---

## Visão geral

O sistema recebe a estrutura de uma molécula (em formato SMILES) e responde a três perguntas em sequência:

1. **A molécula é um inibidor de Cruzipaína ou um composto inativo (decoy)?**
2. **Se for inibidor, qual o seu nível de potência?** (Classe 0 — potente, Classe 1 — moderado, Classe 2 — fraco)
3. **Qual o valor estimado de IC50 em nanomolar?**

Cada pergunta é respondida por um modelo dedicado, e os três operam de forma encadeada.

---

## Dados

- **Fonte:** base de bioatividade ChEMBL, alvo Cruzipaína (CHEMBL3563).
- **Inibidores:** 508 moléculas com IC50 medido, após limpeza e consolidação de duplicatas por mediana.
- **Decoys:** 600 moléculas inativas geradas pela plataforma DUD-E, mantendo proporção próxima de 1:1 com os inibidores.
- **Descritores:** 1875 descritores moleculares (2D e 3D) calculados via PaDEL-Descriptor, reduzidos a **628 atributos** após seleção (remoção de valores ausentes, variância nula e alta correlação).
- **Validação cega:** 3 moléculas (uma de cada classe de potência) foram separadas no início e nunca usadas no treino, reservadas para o teste final do pipeline.

---

## Estrutura do repositório

```
chagas-drug-descovery/
├── app/
│   └── app.py                       # Aplicação Streamlit (pipeline completo)
├── data/
│   ├── raw/                         # Dados brutos do ChEMBL e lotes DUD-E
│   └── processed/                   # Datasets limpos e com descritores
├── models/                          # Modelos e transformadores serializados
├── notebooks/
│   ├── 01_exploratoria.ipynb        # Coleta, limpeza e análise exploratória
│   ├── 02_descritores.ipynb         # Cálculo de descritores e seleção de atributos
│   ├── 03_modelo_classificacao1.ipynb  # Modelo 1 — inibidor vs decoy
│   ├── 04_modelo_classificacao2.ipynb  # Modelo 2 — classificação de potência
│   └── 05_modelo_regressao.ipynb       # Modelo 3 — regressão de IC50
├── requirements.txt
└── README.md
```

---

## Os três modelos

### Modelo 1 — Triagem (inibidor vs decoy)

Classificador binário que separa inibidores de compostos inativos. Foram comparados Logistic Regression, Decision Tree e Random Forest.

O notebook apresenta **duas abordagens**:

- **Referência:** com capping de outliers calculado por grupo, atinge acurácia próxima de 99%.
- **Produção:** o capping por grupo embute informação de classe no pré-processamento e não é replicável para moléculas novas. A versão de produção usa limites de capping fixos (salvos do treino), atinge cerca de 91% de acurácia e é a adotada na aplicação.

Essa distinção corrige um caso de divergência entre treino e produção (*train/serving skew*) e reflete a capacidade real de generalização do modelo.

### Modelo 2 — Classificação de potência (3 classes)

Classifica os inibidores em três níveis de potência, definidos pelos quantis do IC50. Foram comparados Random Forest, HistGradientBoosting, SVM e XGBoost, avaliados por validação cruzada de 10 folds. O melhor desempenho ficou em torno de **76% de acurácia em validação cruzada**, com os modelos baseados em árvore superando o SVM.

### Modelo 3 — Regressão de IC50

Estima o valor contínuo de IC50. A modelagem é feita em **pIC50** (escala logarítmica, padrão em QSAR), o que estabiliza a regressão diante da ampla faixa de valores. Seguindo a metodologia de referência, os atributos são transformados em postos (rank); o modelo prevê o rank do pIC50 e um polinômio converte o resultado de volta, com reconversão final para IC50 em nanomolar.

Desempenho: **R² de aproximadamente 0.92 no espaço pIC50**, com alta precisão nas moléculas potentes, que são as mais relevantes em descoberta de fármacos.

---

## Aplicação interativa

A aplicação em Streamlit reproduz o fluxo completo de predição:

- Entrada por SMILES (até 20 moléculas) ou upload de arquivo CSV/XLSX.
- Cálculo automático de descritores e execução dos três modelos em sequência.
- Resultados com classificação de tipo, classe de potência e IC50 estimado.
- Gráficos de distribuição (tipo e potência) e exportação dos resultados em CSV.
- **Atalho de molécula conhecida:** se o SMILES já existe na base de treino, o resultado é recuperado diretamente, sem reprocessamento pelos modelos.

### Como executar

```
pip install -r requirements.txt
streamlit run app/app.py
```

> O cálculo de descritores via PaDEL requer **Java 8** instalado e acessível no ambiente.

---

## Reprodução do pipeline

Os notebooks devem ser executados na ordem numérica:

1. `01_exploratoria` — coleta, limpeza e separação das moléculas de validação.
2. `02_descritores` — cálculo de descritores PaDEL (com 3D) e seleção de atributos.
3. `03_modelo_classificacao1` — treino e avaliação do Modelo 1.
4. `04_modelo_classificacao2` — treino e avaliação do Modelo 2.
5. `05_modelo_regressao` — treino e avaliação do Modelo 3.

Cada notebook salva seus artefatos em `models/`, consumidos pela aplicação.

---

## Stack técnica

Python, scikit-learn, XGBoost, PaDEL-Descriptor (via padelpy), pandas, NumPy, SciPy, Matplotlib, Seaborn, Streamlit e Plotly.

---

## Limitações e trabalhos futuros

- O conjunto de treino é relativamente pequeno (508 inibidores), o que limita o teto de desempenho e amplia a variância das estimativas; a validação cruzada é usada como métrica de referência por ser mais robusta.
- A reconversão do pIC50 para nanomolar amplifica o erro nas moléculas muito fracas; a avaliação do modelo de regressão é feita no espaço pIC50.
- O reconhecimento de moléculas conhecidas usa comparação textual de SMILES; a canonicalização das estruturas tornaria a busca mais robusta.
- Trabalhos futuros podem incluir ampliação do conjunto de dados, modelagem baseada em grafos moleculares e validação experimental dos candidatos priorizados.

---

## Nota acadêmica

Este projeto foi desenvolvido com fins de estudo e de portfólio, no contexto de modelagem molecular e ciência de dados aplicada à descoberta de fármacos. Reproduz uma metodologia publicada, adaptando-a a um novo alvo (Cruzipaína), e documenta de forma transparente as decisões metodológicas, incluindo as limitações identificadas.
