#!/usr/bin/env python3
"""
mfcc_engine.py
==============
Porta completa dos scripts MFCC_CPCM_A/B/C/D.sh em Python puro.
Sem dependência de bash, awk, sed ou bc — roda em qualquer SO,
inclusive em nuvem (Streamlit Community Cloud).

Lógica de cada modo:
  A → LIG + RES + CAP1 + CAP2          (energia total)
  B → LIG + CAP1 + CAP2                (sem resíduo)
  C → RES + CAP1 + CAP2                (sem ligante)
  D → CAP1 + CAP2                      (só caps)

A energia de interação MFCC é: E_int = E(A) - E(B) - E(C) + E(D)
"""

from __future__ import annotations
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

# ── Tipos ──────────────────────────────────────────────────────────────────────

Mode = Literal["A", "B", "C", "D"]

CHARGED_NEG = {"ASP", "GLU"}
CHARGED_POS = {"LYS", "ARG"}


@dataclass
class Atom:
    element: str
    x: float
    y: float
    z: float


@dataclass
class DatabaseEntry:
    resname: str        # ex: ALA0104
    atom: Atom
    charge: float
    atom_type: str      # ex: CA, N_pep, C_pep, HB1 …


@dataclass
class MFCCContext:
    """Tudo que os scripts precisam para gerar um GJF."""
    database: list[DatabaseEntry]
    residue_list: list[str]          # lista ordenada (lista_residuo.txt)
    ligante: str
    carga_ligante: int
    multiplicidade: int
    bsse: bool
    gaussian_header: str             # conteúdo do gaussian_top.txt
    eps: int


# ── Parser do DATABASE.txt ─────────────────────────────────────────────────────

def parse_database(text: str) -> list[DatabaseEntry]:
    """
    Lê DATABASE.txt e retorna lista de DatabaseEntry.
    Formato de cada linha:
      RESNAME<TAB>ELEM x y z<TAB>CHARGE<TAB>ATOM_TYPE
    """
    entries: list[DatabaseEntry] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        resname   = parts[0].strip()
        coord_str = parts[1].strip()          # "N 37.4780 52.4880 99.5180"
        charge    = float(parts[2].strip())
        atom_type = parts[3].strip()

        coord_parts = coord_str.split()
        if len(coord_parts) < 4:
            continue
        elem = coord_parts[0]
        x, y, z = float(coord_parts[1]), float(coord_parts[2]), float(coord_parts[3])

        entries.append(DatabaseEntry(
            resname=resname,
            atom=Atom(elem, x, y, z),
            charge=charge,
            atom_type=atom_type,
        ))
    return entries


# ── Helpers do DATABASE ────────────────────────────────────────────────────────

def get_pep_c(db: list[DatabaseEntry], resname: str) -> DatabaseEntry | None:
    """Retorna o átomo C_pep do resíduo."""
    for e in db:
        if e.resname == resname and e.atom_type == "C_pep":
            return e
    return None


def get_pep_n(db: list[DatabaseEntry], resname: str) -> DatabaseEntry | None:
    """Retorna o átomo N_pep do resíduo."""
    for e in db:
        if e.resname == resname and e.atom_type == "N_pep":
            return e
    return None


def get_atoms_of(db: list[DatabaseEntry], resname: str) -> list[DatabaseEntry]:
    return [e for e in db if e.resname == resname]


def residue_charge(resname: str, custom: dict[str, int] | None = None) -> int:
    """Determina a carga formal de um resíduo pelos 3 primeiros caracteres."""
    if custom and resname in custom:
        return custom[resname]
    code = resname[:3].upper()
    if code in CHARGED_NEG:
        return -1
    if code in CHARGED_POS:
        return 1
    return 0


# ── Cálculo dos caps de hidrogênio ────────────────────────────────────────────
#
#  Posição do H que substitui a ligação peptídica cortada:
#
#  Cap ANT (C_pep do cap anterior → N_pep do CAP1):
#    H = C_cap1_ant + 0.25 * (N_cap1 - C_cap1_ant)
#
#  Cap POS (N_pep do cap posterior ← C_pep do CAP2):
#    H = N_cap2_pos + 0.20 * (C_cap2 - N_cap2_pos)
#
#  Cap RESANT (C_pep do RESIDUO → N_pep do CAP2):   [scripts B, C, D]
#    H = C_residuo + 0.25 * (N_cap2 - C_residuo)
#
#  Cap RESPOS (N_pep do RESIDUO ← C_pep do CAP1):   [scripts B, C, D]
#    H = N_residuo + 0.25 * (C_cap1 - N_residuo)

def _interp(a: Atom, b: Atom, t: float) -> tuple[float, float, float]:
    """Interpola: a + t*(b - a)"""
    return (
        a.x + t * (b.x - a.x),
        a.y + t * (b.y - a.y),
        a.z + t * (b.z - a.z),
    )


def cap_ant(c_cap1_ant: Atom, n_cap1: Atom) -> tuple[float, float, float]:
    return _interp(c_cap1_ant, n_cap1, 0.25)


def cap_pos(n_cap2_pos: Atom, c_cap2: Atom) -> tuple[float, float, float]:
    return _interp(n_cap2_pos, c_cap2, 0.20)


def cap_resant(c_residuo: Atom, n_cap2: Atom) -> tuple[float, float, float]:
    return _interp(c_residuo, n_cap2, 0.25)


def cap_respos(n_residuo: Atom, c_cap1: Atom) -> tuple[float, float, float]:
    return _interp(n_residuo, c_cap1, 0.25)


# ── Formatação de coordenadas ──────────────────────────────────────────────────

def fmt_atom(elem: str, x: float, y: float, z: float) -> str:
    """Formata linha de átomo no estilo dos scripts: 'H \t x \t y \t z'"""
    return f"{elem} \t {x:.4f} \t {y:.4f} \t {z:.4f}"


def fmt_entry(e: DatabaseEntry) -> str:
    return fmt_atom(e.atom.element, e.atom.x, e.atom.y, e.atom.z)


def fix_leading_dot(text: str) -> str:
    """
    Replica: sed -r 's/(^| |-)(\.[0-9])/\10\2/g'
    Garante que números como -.5 ou .5 virem -0.5 ou 0.5
    """
    return re.sub(r'(^| |-)(\.[0-9])', lambda m: m.group(1) + '0' + m.group(2), text)


# ── Gerador de GJF ────────────────────────────────────────────────────────────

def _neighbors(residue_list: list[str], residue: str):
    """
    Retorna (cap1, cap1_ant, cap2, cap2_pos) para o resíduo.
    cap1     = resíduo anterior
    cap1_ant = anterior ao anterior
    cap2     = resíduo posterior
    cap2_pos = posterior ao posterior
    """
    try:
        idx = residue_list.index(residue)
    except ValueError:
        raise ValueError(f"Resíduo '{residue}' não encontrado em lista_residuo.txt")

    cap1     = residue_list[idx - 1] if idx > 0                          else None
    cap1_ant = residue_list[idx - 2] if idx > 1                          else None
    cap2     = residue_list[idx + 1] if idx < len(residue_list) - 1      else None
    cap2_pos = residue_list[idx + 2] if idx < len(residue_list) - 2      else None

    return cap1, cap1_ant, cap2, cap2_pos


def generate_gjf(
    ctx: MFCCContext,
    residue: str,
    mode: Mode,
    custom_charges: dict[str, int] | None = None,
) -> str:
    """
    Gera o conteúdo de um arquivo GJF para o resíduo e modo dados.
    Replica exatamente a lógica de MFCC_CPCM_A/B/C/D.sh.
    """
    db   = ctx.database
    rlist = ctx.residue_list

    cap1, cap1_ant, cap2, cap2_pos = _neighbors(rlist, residue)

    # ── Coordenadas dos átomos de ligação peptídica ────────────────────────
    def _c(res): return get_pep_c(db, res) if res else None
    def _n(res): return get_pep_n(db, res) if res else None

    C_cap1_ant = _c(cap1_ant)
    N_cap2_pos = _n(cap2_pos)
    C_cap1     = _c(cap1)
    N_cap1     = _n(cap1)
    C_cap2     = _c(cap2)
    N_cap2     = _n(cap2)
    C_residuo  = _c(residue)
    N_residuo  = _n(residue)

    # ── Carga do sistema por modo ──────────────────────────────────────────
    ch_lig  = ctx.carga_ligante
    ch_res  = residue_charge(residue, custom_charges)
    ch_cap1 = residue_charge(cap1, custom_charges) if cap1 else 0
    ch_cap2 = residue_charge(cap2, custom_charges) if cap2 else 0

    charge_map = {
        "A": ch_lig + ch_res + ch_cap1 + ch_cap2,
        "B": ch_lig          + ch_cap1 + ch_cap2,
        "C":          ch_res + ch_cap1 + ch_cap2,
        "D":                   ch_cap1 + ch_cap2,
    }
    carga_sistema = charge_map[mode]

    # ── Descrição do cálculo ───────────────────────────────────────────────
    desc_map = {
        "A": f"Energia LIG+RES+CAP+PC: Lig:{ctx.ligante} + Res:{residue} + Cap1:{cap1} + Cap2:{cap2} ",
        "B": f"Energia LIG+CAP+PC: Lig:{ctx.ligante} + Cap1:{cap1} + Cap2:{cap2} ",
        "C": f"Energia RES+CAP+PC: Res:{residue} + Cap1:{cap1} + Cap2:{cap2}  ",
        "D": f"Energia CAP+PC: Cap1:{cap1} + Cap2:{cap2} ",
    }

    # ── Montar corpo do GJF ────────────────────────────────────────────────
    lines: list[str] = []

    # Cabeçalho Gaussian
    lines.append(ctx.gaussian_header.rstrip())
    lines.append("")
    lines.append(desc_map[mode])
    lines.append("  ")
    lines.append(f"{carga_sistema} {ctx.multiplicidade}")

    # Átomos do ligante
    for e in get_atoms_of(db, ctx.ligante):
        lines.append(fmt_entry(e))

    # Átomos do resíduo (modo A e C)
    if mode in ("A", "C"):
        for e in get_atoms_of(db, residue):
            lines.append(fmt_entry(e))

    # Átomos do CAP1
    for e in get_atoms_of(db, cap1):
        lines.append(fmt_entry(e))

    # Átomos do CAP2
    for e in get_atoms_of(db, cap2):
        lines.append(fmt_entry(e))

    # ── Caps de hidrogênio ─────────────────────────────────────────────────
    # Cap ANT: sempre presente (substitui C_pep do cap1_ant → N do cap1)
    if C_cap1_ant and N_cap1:
        hx, hy, hz = cap_ant(C_cap1_ant.atom, N_cap1.atom)
        lines.append(fmt_atom("H", hx, hy, hz))

    # Cap POS: sempre presente (substitui N_pep do cap2_pos ← C do cap2)
    if N_cap2_pos and C_cap2:
        hx, hy, hz = cap_pos(N_cap2_pos.atom, C_cap2.atom)
        lines.append(fmt_atom("H", hx, hy, hz))

    # Caps extras para o resíduo (modos B, C, D — onde o resíduo foi removido)
    if mode in ("B", "C", "D"):
        if C_residuo and N_cap2:
            hx, hy, hz = cap_resant(C_residuo.atom, N_cap2.atom)
            lines.append(fmt_atom("H", hx, hy, hz))
        if N_residuo and C_cap1:
            hx, hy, hz = cap_respos(N_residuo.atom, C_cap1.atom)
            lines.append(fmt_atom("H", hx, hy, hz))

    lines.append("  ")
    lines.append(f"eps={ctx.eps}")
    lines.append("  ")

    raw = "\n".join(lines)
    raw = fix_leading_dot(raw)

    # Remover linhas com coordenadas 900.0000 (gaps artificiais)
    raw = "\n".join(
        l for l in raw.splitlines()
        if "900.0000" not in l
    )

    return raw


# ── API pública ────────────────────────────────────────────────────────────────

def build_gaussian_header(
    mem: str, nproc: int, method: str, basis_set: str,
    scrf: str = "cpcm,solvent=water,read",
    extra: str = "density=current NoSymmetry scf=(xqc,maxcycles=1024)",
) -> str:
    return (
        f"%mem={mem}\n"
        f"%nproc={nproc}\n"
        f"#p {method}/{basis_set} scrf=({scrf}) {extra}"
    )


def run_mfcc_python(
    database_text: str,
    residue_list_text: str,
    mfcc_residuos_text: str,
    params: dict,
    log_fn=print,
) -> dict[str, str]:
    """
    Executa o pipeline MFCC completo em Python puro.

    Parâmetros
    ----------
    database_text      : conteúdo de DATABASE.txt
    residue_list_text  : conteúdo de lista_residuo.txt
    mfcc_residuos_text : conteúdo de MFCC_RESIDUOS.txt
    params             : dicionário com as configurações do usuário
    log_fn             : função de log (padrão: print)

    Retorna
    -------
    dict: {nome_arquivo: conteúdo_gjf}
    """
    db      = parse_database(database_text)
    rlist   = [r.strip() for r in residue_list_text.splitlines() if r.strip()]
    targets = [r.strip() for r in mfcc_residuos_text.splitlines() if r.strip()]

    header = build_gaussian_header(
        mem=params["mem"],
        nproc=params["nproc"],
        method=params["method"],
        basis_set=params["basis_set"],
    )

    ctx = MFCCContext(
        database=db,
        residue_list=rlist,
        ligante=params["ligand_resname"],
        carga_ligante=params["carga_ligante"],
        multiplicidade=params["multiplicidade"],
        bsse=params.get("bsse", False),
        gaussian_header=header,
        eps=params["eps"],
    )

    custom_charges = params.get("charged_residues", {})
    results: dict[str, str] = {}
    errors: list[str] = []

    for i, residue in enumerate(targets, 1):
        for mode in ("A", "B", "C", "D"):
            fname = f"{ctx.ligante}_{residue}_{mode}.gjf"
            try:
                content = generate_gjf(ctx, residue, mode, custom_charges)
                results[fname] = content
            except Exception as e:
                errors.append(f"{fname}: {e}")
                log_fn(f"    [FALHA] {fname}: {e}")

        log_fn(f"  [{i:3d}/{len(targets)}] OK  {residue}")

    if errors:
        log_fn(f"\n[AVISO] {len(errors)} erro(s) encontrado(s)")

    log_fn(f"\n[OK] {len(results)} GJFs gerados")
    return results
