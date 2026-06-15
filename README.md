# Chagas Drug Discovery - Triagem Virtual de Inibidores de Cruzipain

Pipeline de aprendizado de maquina para a descoberta de inibidores da **Cruzipain** (Cruzain), a principal cisteino-protease do *Trypanosoma cruzi*, agente causador da doenca de Chagas. O projeto reproduz e adapta uma metodologia de triagem virtual em tres etapas, originalmente desenvolvida para um alvo de Alzheimer, aplicando-a a um novo alvo terapeutico.

A doenca de Chagas afeta milhoes de pessoas, sobretudo na America Latina, e dispoe de poucas opcoes terapeuticas. A Cruzipain e essencial para a sobrevivencia e a replicacao do parasita, o que a torna um alvo molecular promissor.

---

## Visao geral

O sistema recebe a estrutura de uma molecula (em formato SMILES) e responde a tres perguntas em sequencia:

1. **A molecula e um inibidor de Cruzipain ou um composto inativo (decoy)?**
2. **Se for inibidor, qual o seu nivel de potencia?** (Classe 0 - potente, Classe 1 - moderado, Classe 2 - fraco)
3. **Qual o valor estimado de IC50 em nanomolar?**

Cada pergunta e respondida por um modelo dedicado, e os tres operam de forma encadeada.

---

## Dados

- **Fonte:** base de bioatividade ChEMBL, alvo Cruzipain (CHEMBL3563).
- **Inibidores:** 508 moleculas com IC50 medido, apos limpeza e consolidacao de duplicatas por mediana.
- **Decoys:** 600 moleculas inativas geradas pela plataforma DUD-E, mantendo proporcao proxima de 1:1 com os inibidores.
- **Descritores:** 1875 descritores moleculares (2D e 3D) calculados via PaDEL-Descriptor, reduzidos a **628 atributos** apos selecao (remocao de valores ausentes, variancia nula e alta correlacao).
- **Validacao cega:** 3 moleculas (uma de cada classe de potencia) foram separadas no inicio e nunca usadas no treino, reservadas para o teste final do pipeline.

---

## Estrutura do repositorio

```
chagas-drug-descovery/
├── app/
│   └── app.py                       # Aplicacao Streamlit (pipeline completo)
├── data/
│   ├── raw/                         # Dados brutos do ChEMBL e lotes DUD-E
│   └── processed/                   # Datasets limpos e com descritores
├── models/                          # Modelos e transformadores serializados
├── notebooks/
│   ├── 01_exploratoria.ipynb        # Coleta, limpeza e analise exploratoria
│   ├── 02_descritores.ipynb         # Calculo de descritores e selecao de atributos
│   ├── 03_modelo_classificacao1.ipynb  # Modelo 1 - inibidor vs decoy
│   ├── 04_modelo_classificacao2.ipynb  # Modelo 2 - classificacao de potencia
│   └── 05_modelo_regressao.ipynb       # Modelo 3 - regressao de IC50
├── requirements.txt
└── README.md
```

---

## Os tres modelos

### Modelo 1 - Triagem (inibidor vs decoy)

Classificador binario que separa inibidores de compostos inativos. Foram comparados Logistic Regression, Decision Tree e Random Forest.

O notebook apresenta **duas abordagens**:
- **Referencia:** com capping de outliers calculado por grupo, atinge acuracia proxima de 99%.
- **Producao:** o capping por grupo embute informacao de classe no pre-processamento e nao e replicavel para moleculas novas. A versao de producao usa limites de capping fixos (salvos do treino), atinge cerca de 91% de acuracia e e a adotada na aplicacao.

Essa distincao corrige um caso de divergencia entre treino e producao (*train/serving skew*) e reflete a capacidade real de generalizacao do modelo.

### Modelo 2 - Classificacao de potencia (3 classes)

Classifica os inibidores em tres niveis de potencia, definidos pelos quantis do IC50. Foram comparados Random Forest, HistGradientBoosting, SVM e XGBoost, avaliados por validacao cruzada de 10 folds. O melhor desempenho ficou em torno de **76% de acuracia em validacao cruzada**, com os modelos baseados em arvore superando o SVM.

### Modelo 3 - Regressao de IC50

Estima o valor continuo de IC50. A modelagem e feita em **pIC50** (escala logaritmica, padrao em QSAR), o que estabiliza a regressao diante da ampla faixa de valores. Seguindo a metodologia de referencia, os atributos sao transformados em postos (rank); o modelo preve o rank do pIC50 e um polinomio converte o resultado de volta, com reconversao final para IC50 em nanomolar.

Desempenho: **R2 de aproximadamente 0.92 no espaco pIC50**, com alta precisao nas moleculas potentes, que sao as mais relevantes em descoberta de farmacos.

---

## Aplicacao interativa

A aplicacao em Streamlit reproduz o fluxo completo de predicao:

- Entrada por SMILES (ate 20 moleculas) ou upload de arquivo CSV/XLSX.
- Calculo automatico de descritores e execucao dos tres modelos em sequencia.
- Resultados com classificacao de tipo, classe de potencia e IC50 estimado.
- Graficos de distribuicao (tipo e potencia) e exportacao dos resultados em CSV.
- **Atalho de molecula conhecida:** se o SMILES ja existe na base de treino, o resultado e recuperado diretamente, sem reprocessamento pelos modelos.

### Como executar

```bash
pip install -r requirements.txt
streamlit run app/app.py
```

> O calculo de descritores via PaDEL requer **Java 8** instalado e acessivel no ambiente.

---

## Reproducao do pipeline

Os notebooks devem ser executados na ordem numerica:

1. `01_exploratoria` - coleta, limpeza e separacao das moleculas de validacao.
2. `02_descritores` - calculo de descritores PaDEL (com 3D) e selecao de atributos.
3. `03_modelo_classificacao1` - treino e avaliacao do Modelo 1.
4. `04_modelo_classificacao2` - treino e avaliacao do Modelo 2.
5. `05_modelo_regressao` - treino e avaliacao do Modelo 3.

Cada notebook salva seus artefatos em `models/`, consumidos pela aplicacao.

---

## Stack tecnica

Python, scikit-learn, XGBoost, PaDEL-Descriptor (via padelpy), pandas, NumPy, SciPy, Matplotlib, Seaborn, Streamlit e Plotly.

---

## Limitacoes e trabalhos futuros

- O conjunto de treino e relativamente pequeno (508 inibidores), o que limita o teto de desempenho e amplia a variancia das estimativas; a validacao cruzada e usada como metrica de referencia por ser mais robusta.
- A reconversao do pIC50 para nanomolar amplifica o erro nas moleculas muito fracas; a avaliacao do modelo de regressao e feita no espaco pIC50.
- O reconhecimento de moleculas conhecidas usa comparacao textual de SMILES; a canonicalizacao das estruturas tornaria a busca mais robusta.
- Trabalhos futuros podem incluir ampliacao do conjunto de dados, modelagem baseada em grafos moleculares e validacao experimental dos candidatos priorizados.

---

## Nota academica

Este projeto foi desenvolvido com fins de estudo e de portfolio, no contexto de modelagem molecular e ciencia de dados aplicada a descoberta de farmacos. Reproduz uma metodologia publicada, adaptando-a a um novo alvo (Cruzipain), e documenta de forma transparente as decisoes metodologicas, incluindo as limitacoes identificadas.
