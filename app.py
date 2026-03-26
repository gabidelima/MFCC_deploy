#!/usr/bin/env python3
"""
app.py — MFCC Pipeline · Interface Streamlit
100% Python puro — sem bash, awk, sed ou bc.
Deploy direto no Streamlit Community Cloud.
Execute com: streamlit run app.py
"""

import io
import tempfile
import traceback
import zipfile
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="MFCC Pipeline", page_icon="⚗️",
                   layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@400;600;800&display=swap');
html, body, [class*="css"] { font-family: 'Syne', sans-serif; }
.stApp { background: #0d1117; color: #c9d1d9; }
[data-testid="stSidebar"] { background: #0d1117; border-right: 1px solid #21262d; }
.main-title { font-family: 'Syne', sans-serif; font-weight: 800; font-size: 2.2rem;
    color: #f0f6fc; letter-spacing: -0.03em; line-height: 1.1; margin-bottom: 0.2rem; }
.main-sub { font-family: 'JetBrains Mono', monospace; font-size: 0.72rem; color: #388bfd;
    letter-spacing: 0.15em; text-transform: uppercase; margin-bottom: 1.5rem; }
.stTextInput input, .stNumberInput input {
    background: #0d1117 !important; border: 1px solid #30363d !important;
    border-radius: 6px !important; color: #c9d1d9 !important;
    font-family: 'JetBrains Mono', monospace !important; }
.stSelectbox > div > div { background: #0d1117 !important; border: 1px solid #30363d !important; color: #c9d1d9 !important; }
.stButton > button { background: #238636 !important; color: #fff !important;
    border: 1px solid #2ea043 !important; border-radius: 6px !important;
    font-family: 'Syne', sans-serif !important; font-weight: 600 !important;
    width: 100%; transition: background 0.15s !important; }
.stButton > button:hover { background: #2ea043 !important; }
.log-box { background: #010409; border: 1px solid #21262d; border-radius: 6px;
    padding: 1rem 1.25rem; font-family: 'JetBrains Mono', monospace;
    font-size: 0.72rem; color: #7ee787; line-height: 1.8;
    max-height: 360px; overflow-y: auto; white-space: pre-wrap; }
hr { border-color: #21262d; }
label { color: #8b949e !important; font-size: 0.82rem !important; }
</style>
""", unsafe_allow_html=True)

METHODS = ["B97D","B3LYP","M062X","wB97XD","PBE0","CAM-B3LYP","B3PW91","TPSS","BP86","HF"]
BASIS_SETS = ["6-311+G(d,p)","6-31G(d)","6-31+G(d,p)","6-311G(d,p)","6-311++G(d,p)",
              "cc-pVDZ","cc-pVTZ","aug-cc-pVDZ","def2-SVP","def2-TZVP"]

def build_zip(gjf_dict, db_text=None, csv_text=None, ligand_resname="", ligand_resid=""):
    """
    Cria ZIP com estrutura:
    - LIGANTE_NUMERO/ (pasta com todos os .gjf)
    - DATABASE.txt (na raiz)
    - residuos_raio.csv (na raiz)
    
    Exemplo: GW11197/
    """
    # Criar nome da pasta: LIGANTE + RESID (ex: GW11197)
    if ligand_resid:
        folder_name = f"{ligand_resname}{ligand_resid}"
    else:
        folder_name = "gjf"
    
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        # Adicionar GJFs em pasta nomeada
        for fname, content in sorted(gjf_dict.items()):
            zf.writestr(f"{folder_name}/{fname}", content)
        
        # Adicionar DATABASE.txt na raiz
        if db_text:
            zf.writestr("DATABASE.txt", db_text)
        
        # Adicionar residuos_raio.csv na raiz
        if csv_text:
            zf.writestr("residuos_raio.csv", csv_text)
    
    return buf.getvalue()

def pdb_stats(pdb_text):
    lines = pdb_text.splitlines()
    atoms   = [l for l in lines if l.startswith("ATOM")]
    hetatms = [l for l in lines if l.startswith("HETATM")]
    res_set = {l[17:20].strip() for l in atoms}
    return len(atoms), len(hetatms), len(res_set)

def ligand_in_pdb(pdb_text, ligand):
    return any(l.startswith("HETATM") and l[17:20].strip() == ligand
               for l in pdb_text.splitlines())

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown('<div class="main-title">⚗️ MFCC</div>', unsafe_allow_html=True)
    st.markdown('<div class="main-sub">Pipeline Automation</div>', unsafe_allow_html=True)
    st.markdown("---")
    st.markdown("#### 🧬 Ligante")
    ligand_resname = st.text_input("Código 3 letras do ligante", value="AGO", max_chars=10,
        help="Nome do resíduo HETATM no PDB").strip().upper()
    col1, col2 = st.columns(2)
    with col1:
        chain = st.text_input("Cadeia", value="", max_chars=2,
            help="Deixe VAZIO para aceitar todas as cadeias. A proteína pode estar em cadeia diferente do ligante!").strip().upper() or None
    with col2:
        carga_ligante = st.number_input("Carga", value=0, min_value=-10, max_value=10, step=1)
    multiplicidade = st.number_input("Multiplicidade", value=1, min_value=1, max_value=5, step=1)
    st.markdown("---")
    st.markdown("#### ✂️ Fragmentação")
    cutoff = st.number_input("Raio de corte (Å)", value=8.0,
        min_value=1.0, max_value=20.0, step=0.5, format="%.1f")
    st.markdown("---")
    st.markdown("#### ⚙️ Gaussian")
    method    = st.selectbox("Método (DFT/HF)", METHODS, index=0)
    basis_set = st.selectbox("Basis set", BASIS_SETS, index=0)
    col3, col4 = st.columns(2)
    with col3:
        mem = st.text_input("Memória", value="16GB")
    with col4:
        nproc = st.number_input("nproc", value=8, min_value=1, max_value=512, step=1)
    eps  = st.number_input("Constante dielétrica (ε)", value=40, min_value=1, max_value=200, step=1)
    bsse = st.toggle("Correção BSSE", value=False)
    st.markdown("---")
    st.markdown('<span style="font-family:\'JetBrains Mono\',monospace;font-size:0.62rem;color:#484f58;">'
                'MFCC Automation Pipeline<br>github.com/seu-usuario/mfcc-automation</span>',
                unsafe_allow_html=True)

# ── Área principal ────────────────────────────────────────────────────────────
st.markdown('<div class="main-title">MFCC Automation Pipeline</div>', unsafe_allow_html=True)
st.markdown('<div class="main-sub">Molecular Fractionation with Conjugate Caps — PDB → GJF</div>',
            unsafe_allow_html=True)

c1,c2,c3,c4,c5 = st.columns(5)
c1.metric("Ligante",  ligand_resname)
c2.metric("Cutoff",   f"{cutoff} Å")
c3.metric("Método",   method)
c4.metric("Basis",    basis_set)
c5.metric("ε",        eps)
st.markdown("---")

# ── Step 1: Upload ───────────────────────────────────────────────────────────
st.markdown("### 1 — Upload do PDB")
uploaded = st.file_uploader("Arraste ou selecione o PDB preparado (com hidrogênios)", type=["pdb"])

if not uploaded:
    st.info("⬆️ Faça o upload de um arquivo PDB para começar.")
    st.markdown("""
    <div style="background:#161b22;border:1px dashed #30363d;border-radius:8px;
                padding:3rem;text-align:center;margin-top:1rem;">
        <div style="font-size:3rem;margin-bottom:1rem;">⚗️</div>
        <div style="font-family:'JetBrains Mono',monospace;font-size:0.75rem;color:#484f58;line-height:2.2;">
            1. Configure os parâmetros na barra lateral<br>
            2. Faça upload do PDB preparado (com H)<br>
            3. Clique em Executar Pipeline<br>
            4. Baixe os GJFs prontos para o Gaussian
        </div>
    </div>""", unsafe_allow_html=True)
    st.stop()

pdb_bytes = uploaded.read()
pdb_text  = pdb_bytes.decode("utf-8", errors="ignore")
n_atoms, n_hetatm, n_res = pdb_stats(pdb_text)
ca, cb, cc = st.columns(3)
ca.metric("Átomos ATOM",     n_atoms)
cb.metric("Átomos HETATM",   n_hetatm)
cc.metric("Resíduos únicos", n_res)

if ligand_in_pdb(pdb_text, ligand_resname):
    st.success(f"✓ Ligante **{ligand_resname}** encontrado no PDB")
else:
    st.warning(f"⚠️ Ligante **{ligand_resname}** não encontrado nos HETATM. Verifique o código.")

st.markdown("---")
st.markdown("### 2 — Executar Pipeline")

with st.expander("📋 Revisar configurações", expanded=False):
    st.code(f"""ligand_resname : {ligand_resname}
chain          : {chain or 'todas'}
carga_ligante  : {carga_ligante}
multiplicidade : {multiplicidade}
cutoff         : {cutoff} Å
method         : {method}
basis_set      : {basis_set}
mem            : {mem}
nproc          : {nproc}
eps            : {eps}
bsse           : {bsse}""", language="yaml")

if st.button("▶ Executar MFCC Pipeline", type="primary"):
    st.markdown("### 3 — Execução")
    log_box  = st.empty()
    progress = st.progress(0, text="Iniciando…")
    log_lines = []

    def log(msg):
        log_lines.append(msg)
        log_box.markdown(
            f'<div class="log-box">{"<br>".join(log_lines)}</div>',
            unsafe_allow_html=True)

    try:
        with tempfile.TemporaryDirectory(prefix="mfcc_") as tmp_str:
            tmp = Path(tmp_str)
            pdb_path = tmp / "input.pdb"
            pdb_path.write_bytes(pdb_bytes)

            # Etapa 1 — DATABASE
            log("[1/3] Gerando DATABASE.txt e lista_residuo.txt …")
            progress.progress(5, text="Gerando DATABASE…")
            from generate_database import parse_pdb, write_outputs as write_db
            db_lines, res_list = parse_pdb(str(pdb_path), chain, ligand_resname)
            if not res_list:
                st.error("Nenhum resíduo proteico encontrado. Verifique cadeia e PDB.")
                st.stop()
            db_path, lr_path = write_db(db_lines, res_list, tmp)
            db_text = db_path.read_text()
            rl_text = lr_path.read_text()
            n_lig  = sum(1 for r, _ in db_lines if r is None)
            n_prot = sum(1 for r, _ in db_lines if r is not None)
            log(f"    {n_lig} átomos ligante | {n_prot} átomos proteicos | {len(res_list)} resíduos")

            # Etapa 2 — Fragmentação
            # Etapa 2 — Fragmentação
            log(f"[2/3] Selecionando resíduos (cutoff={cutoff} Å) …")
            progress.progress(20, text="Calculando distâncias…")
            from fragment_by_distance import find_residues_in_cutoff, write_outputs as write_frag
            try:
                results_dist, ligand_resid = find_residues_in_cutoff(str(pdb_path), ligand_resname, chain, float(cutoff))
            except ValueError as e:
                log(f"\n❌ ERRO: {str(e)}")
                st.error(str(e))
                st.stop()
            
            if not results_dist:
                error_msg = f"❌ Erro: Nenhum resíduo dentro de {cutoff} Å encontrado.\n\n"
                error_msg += "Possíveis causas:\n"
                error_msg += "1. ⚠️ **Ligante e proteína em cadeias diferentes**\n"
                error_msg += "   → Deixe a 'Cadeia' VAZIA para aceitar todas\n"
                error_msg += "2. Cutoff muito pequeno para a estrutura\n"
                error_msg += "3. Código 3-letras do ligante incorreto\n\n"
                error_msg += "💡 Dica: Se você selecionou uma cadeia específica,\n"
                error_msg += "    **deixe em branco** — a proteína pode estar em\n"
                error_msg += "    cadeia diferente do ligante!"
                st.error(error_msg)
                st.stop()
            
            csv_path, mfcc_path = write_frag(results_dist, tmp)
            csv_text = csv_path.read_text()
            mr_text = mfcc_path.read_text()
            targets = [r.strip() for r in mr_text.splitlines() if r.strip()]
            log(f"    {len(targets)} resíduos selecionados")

            # Etapa 3 — GJFs
            log(f"[3/3] Gerando GJFs (4 modos × {len(targets)} = {4*len(targets)} arquivos) …")
            progress.progress(30, text="Gerando GJFs…")

            params = {
                "ligand_resname": ligand_resname, "chain": chain,
                "cutoff": cutoff, "method": method, "basis_set": basis_set,
                "mem": mem, "nproc": int(nproc), "eps": int(eps),
                "bsse": bsse, "carga_ligante": int(carga_ligante),
                "multiplicidade": int(multiplicidade),
            }

            import re as _re
            def log_progress(msg):
                m = _re.search(r'\[\s*(\d+)/\s*(\d+)\]', msg)
                if m:
                    done, total = int(m.group(1)), int(m.group(2))
                    pct = 30 + int(done / total * 65)
                    progress.progress(pct, text=f"Resíduo {done}/{total}…")
                log(msg)

            from mfcc_engine import run_mfcc_python
            gjf_dict = run_mfcc_python(db_text, rl_text, mr_text, params, log_fn=log_progress)
            progress.progress(100, text="Concluído!")
            log(f"\n{'='*48}")
            log(f"[CONCLUÍDO] {len(gjf_dict)} GJFs gerados com sucesso")

        # Step 4 — Download
        st.markdown("---")
        st.markdown("### 4 — Download")
        st.success(f"✅ **{len(gjf_dict)} arquivos GJF** gerados com sucesso!")

        sample = sorted(gjf_dict.keys())[0]
        with st.expander(f"👁 Preview: {sample}"):
            st.code(gjf_dict[sample], language="text")

        st.download_button(
            label=f"⬇ Baixar todos os GJFs ({len(gjf_dict)} arquivos) — ZIP",
            data=build_zip(gjf_dict, db_text, csv_text, ligand_resname, ligand_resid),
            file_name=f"{ligand_resname}_MFCC_gjf.zip",
            mime="application/zip",
        )

        with st.expander("📊 Resíduos selecionados e distâncias ao ligante"):
            import pandas as pd
            df = pd.DataFrame(results_dist, columns=["Distância (Å)", "Resíduo"])
            df["Distância (Å)"] = df["Distância (Å)"].astype(float).round(2)
            st.dataframe(df, use_container_width=True, hide_index=True)

    except Exception as e:
        st.error(f"Erro inesperado: {e}")
        st.code(traceback.format_exc())
