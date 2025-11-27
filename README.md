# 🧬 LDPC_PCOM — Simulações e Implementações de Códigos LDPC

![Status](https://img.shields.io/badge/Status-Ativo-blueviolet?style=for-the-badge)
![MATLAB](https://img.shields.io/badge/MATLAB-OK-blue?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-Implementação%20Ativa-yellow?style=for-the-badge)
![Projeto](https://img.shields.io/badge/Tipo-Acadêmico%20%2F%20Pesquisa-green?style=for-the-badge)

Repositório dedicado ao estudo, implementação e simulação de **Códigos LDPC (Low-Density Parity-Check)** para a disciplina **PCOM — Processamento e Comunicação de Sinais**.

Este projeto inclui implementações em MATLAB e Python do algoritmo **Belief Propagation (SPA — Sum-Product Algorithm)** aplicados ao canal AWGN com mapeamento BPSK.

---

## 🎯 Objetivos do Projeto

- Implementar um código LDPC básico e validá-lo em canal ruidoso.
- Comparar **MATLAB x Python** na simulação do Sum-Product Algorithm.
- Avaliar BER × SNR para diferentes condições de canal.
- Criar uma base modular e extensível para:
  - QC-LDPC
  - Min-Sum Algorithm
  - Paralelização
  - Matrizes H de padrões reais (802.11n, WiMAX, DVB-S2)

---

## 📁 Estrutura do Repositório

```plaintext
LDPC_PCOM/
│
├── matlab/               # Versão MATLAB/Octave do simulador LDPC
│   └── ldpc_demo.m
│
├── python/               # Versão Python do simulador LDPC
│   ├── ldpc_bp.py
│   └── requirements.txt
│
├── docs/                 # Materiais complementares (PDFs, referências)
│
├── media/                # Figuras, diagramas e gráficos gerados
│
├── tests/                # Scripts de teste (unitários ou funcionais)
│
├── results/              # Saída de experimentos, logs, curvas BER etc.
│
├── GUIA.md               # Guia de commits e boas práticas
├── NOTES.md              # Anotações gerais
├── TESTING.md            # Como validar cada módulo
└── README.md             # Este arquivo
```

---

## 🧪 Tecnologias e Metodologia

### 🔹 Canal e Modulação

- **Modulação**: BPSK
- **Canal**: AWGN
- Cálculo direto dos LLRs do canal usando:

$$
L_{ch}(v) = \frac{2 y_v}{\sigma^2}
$$

- Relação entre \(E_b/N_0\) e a variância do ruído configurada conforme o código.

---

### 🔹 Estrutura do Código LDPC

- **Matriz geradora**:

$$
G = [I_k \;\; P]
$$

- **Matriz de verificação de paridade**:

$$
H = [P^T \;\; I_r]
$$

- Construção manual esparsa para experimentos didáticos.
- Representação computacional:
  - **Lista de adjacência** para o Grafo de Tanner.
  - Compatível com MATLAB e Python.
  - Perfis de nós de checagem (CN) e nós variáveis (VN).

---

### 🔹 Decodificação — SPA (Sum-Product Algorithm)

Implementação completa do algoritmo iterativo de passagem de mensagens:

1. **Mensagens VN → CN**

   - Cada variável envia seus LLRs para os nós de checagem vizinhos.

2. **Mensagens CN → VN**

$$
r_{m,v} =
2 \cdot \operatorname{atanh}\left(
\prod_{v' \neq v} \tanh\left(\frac{q_{m,v'}}{2}\right)
\right)
$$

3. **Combinação dos LLRs**

$$
L_{ap}(v) =
L_{ch}(v) + \sum_{m \in \mathcal{N}(v)} r_{m,v}
$$

4. **Decisão dura**

- Bit = 0 se $L_{ap}(v) > 0$, caso contrário 1.

5. **Convergência**

- Se $H \hat{c}^T = 0$, a decodificação finaliza antes do limite de iterações.

---

## ▶️ Execução das Simulações

### MATLAB

```bash
cd matlab
run ldpc_demo.m
```

### Python

```bash
cd python
pip install -r requirements.txt
python ldpc_bp.py
```

---
