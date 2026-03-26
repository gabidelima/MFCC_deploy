#!/usr/bin/env python3
"""
fragment_by_distance.py
=======================
Seleciona resíduos proteicos dentro de um raio de corte ao redor do ligante
e gera mfcc_residuos.txt com os nomes formatados (4 dígitos).

Uso:
    python fragment_by_distance.py --pdb proteina.pdb --ligand AGO --cutoff 8.0

Saída:
    residuos_raio.csv   — tabela distância × resíduo (para conferência)
    mfcc_residuos.txt   — resíduos selecionados, nomes com 4 dígitos

Dependências:
    pip install biopython
"""

import csv
import sys
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


# ---------------------------------------------------------------------------
# Funções auxiliares
# ---------------------------------------------------------------------------

def format_resname_4digits(resname: str, resid: int) -> str:
    """Formata RES + número com exatamente 4 dígitos: ALA0027, GLU0104."""
    return f"{resname.strip()}{resid:04d}"


def min_distance_residue_ligand(residue, ligand) -> float:
    """
    Calcula a menor distância átomo-a-átomo entre um resíduo proteico
    e o ligante.
    """
    min_dist = None
    for lig_atom in ligand.get_atoms():
        for prot_atom in residue.get_atoms():
            d = prot_atom - lig_atom
            if min_dist is None or d < min_dist:
                min_dist = d
    return min_dist


# ---------------------------------------------------------------------------
# Lógica principal
# ---------------------------------------------------------------------------

def find_residues_in_cutoff(pdb_file: str, ligand_resname: str,
                             chain_id=None, cutoff: float = 8.0):
    """
    Percorre o PDB e retorna lista de (distância, label_formatado)
    para todos os resíduos dentro do raio de corte.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mol", pdb_file)
    model = structure[0]

    # --- Encontrar o ligante ---
    ligand = None
    all_hetatm = []
    
    for residue in model.get_residues():
        het, _, _ = residue.get_id()
        resname = residue.get_resname().strip()
        if het.strip():  # É um HETATM
            all_hetatm.append(resname)
            if resname == ligand_resname:
                ligand = residue
                break

    if ligand is None:
        error_msg = f"⚠️ Ligante '{ligand_resname}' não encontrado no PDB.\n"
        if all_hetatm:
            unique_hetatm = sorted(set(all_hetatm))
            error_msg += f"   HETATM disponíveis: {', '.join(unique_hetatm)}\n"
        error_msg += f"   Verifique o nome do resíduo (campo HETATM, coluna 18-20)."
        raise ValueError(error_msg)

    print(f"[INFO] Ligante encontrado: {ligand.get_resname()} "
          f"(chain {ligand.get_parent().get_id()}, "
          f"resid {ligand.get_id()[1]})")

    # --- Calcular distâncias ---
    results = []

    for chain in model.get_chains():
        if chain_id and chain.get_id() != chain_id:
            continue
        for residue in chain.get_residues():
            het, resid, _ = residue.get_id()
            if het.strip():   # pular HETATM (incluindo o próprio ligante)
                continue
            if not is_aa(residue, standard=False):
                continue

            dist = min_distance_residue_ligand(residue, ligand)
            if dist is not None and dist <= cutoff:
                label = format_resname_4digits(residue.get_resname(), resid)
                results.append((dist, label))

    results.sort(key=lambda x: x[0])
    return results


def write_outputs(results, output_dir: Path):
    """Escreve residuos_raio.csv e mfcc_residuos.txt."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path    = output_dir / "residuos_raio.csv"
    mfcc_path   = output_dir / "mfcc_residuos.txt"

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Distância (Å)", "Resíduo"])
        for dist, label in results:
            writer.writerow([f"{dist:.4f}", label])

    with open(mfcc_path, "w") as f:
        for _, label in results:
            f.write(label + "\n")

    return csv_path, mfcc_path
