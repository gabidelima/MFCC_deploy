#!/usr/bin/env python3
"""
fragment_by_distance.py
=======================
Encontra resíduos dentro de um raio de corte (cutoff) a partir do ligante.
Usa a distância euclidiana entre centros de massa ou átomos representativos.
"""

from pathlib import Path
from typing import Optional
import math


def parse_pdb_ligand(pdb_file: str, ligand_resname: str, chain: Optional[str] = None) -> dict:
    """
    Extrai coordenadas de todos os átomos do ligante.
    Retorna: {"ligand_atoms": [(x, y, z), ...], "residues": {res_id: [(x,y,z), ...]}}
    """
    pdb_lines = Path(pdb_file).read_text().splitlines()
    
    ligand_atoms = []
    residues = {}
    all_hetatm = set()  # Para debug
    
    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        
        record = line[0:6].strip()
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain_id = line[21].strip()
        
        # Coletar todos os HETATM para debug
        if record == "HETATM":
            all_hetatm.add(resname)
        
        try:
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except (ValueError, IndexError):
            continue
        
        # Filtrar por cadeia (filtro suave: se chain não é fornecido, aceita tudo)
        if chain and chain_id != chain:
            continue
        
        # Ligante: procura QUALQUER HETATM (não apenas o especificado)
        # Se nenhum HETATM for encontrado com o nome exato, lista todas as opções
        if record == "HETATM":
            if resname == ligand_resname:
                ligand_atoms.append((x, y, z))
        
        # Proteína (ATOM records)
        elif record == "ATOM":
            res_id = f"{resname}{res_seq}"
            if res_id not in residues:
                residues[res_id] = []
            residues[res_id].append((x, y, z))
    
    return {
        "ligand_atoms": ligand_atoms,
        "residues": residues,
        "all_hetatm": list(all_hetatm),  # Para debug
    }


def distance_3d(p1: tuple, p2: tuple) -> float:
    """Distância euclidiana entre dois pontos 3D."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)


def find_residues_in_cutoff(
    pdb_file: str,
    ligand_resname: str,
    chain: Optional[str] = None,
    cutoff: float = 8.0,
    log_fn = None  # Função de logging opcional
) -> list:
    """
    Encontra resíduos dentro de cutoff Ångströms do ligante.
    Usa a DISTÂNCIA MÍNIMA entre qualquer átomo do resíduo e qualquer átomo do ligante.
    
    Retorna:
        lista de tuplas (distância, resíduo_id), ordenada por distância
    """
    data = parse_pdb_ligand(pdb_file, ligand_resname, chain)
    
    # Debug: verificar se o ligante foi encontrado
    if not data["ligand_atoms"]:
        msg = f"\n⚠️ AVISO: Nenhum átomo do ligante '{ligand_resname}' encontrado!"
        if data["all_hetatm"]:
            msg += f"\n   HETATM disponíveis no PDB: {', '.join(sorted(data['all_hetatm']))}"
            msg += f"\n   Verifique se '{ligand_resname}' está correto."
        if log_fn:
            log_fn(msg)
        return []
    
    if log_fn:
        log_fn(f"[DEBUG] Ligante encontrado: {len(data['ligand_atoms'])} átomos")
        log_fn(f"[DEBUG] Proteína: {len(data['residues'])} resíduos únicos")
    
    results = []
    
    # Para cada resíduo, calcula a distância mínima até o ligante
    for res_id, res_coords in data["residues"].items():
        min_distance = float("inf")
        
        # Testa todos os pares átomo-ligante
        for atom_lig in data["ligand_atoms"]:
            for atom_res in res_coords:
                dist = distance_3d(atom_lig, atom_res)
                if dist < min_distance:
                    min_distance = dist
                
                # Early exit se já passou o cutoff
                if min_distance > cutoff:
                    break
            if min_distance > cutoff:
                break
        
        if min_distance <= cutoff:
            results.append((min_distance, res_id))
    
    # Ordenar por distância
    results.sort(key=lambda x: x[0])
    
    if log_fn and results:
        log_fn(f"[DEBUG] Encontrados {len(results)} resíduos dentro de {cutoff} Å")
    
    return results


def write_outputs(results_dist: list, output_dir: Path) -> tuple:
    """
    Escreve residuos_raio.csv e mfcc_residuos.txt.
    
    Retorna:
        (csv_path, residues_path)
    """
    output_dir = Path(output_dir)
    
    # residuos_raio.csv: distância,resíduo
    csv_lines = ["Distância (Å),Resíduo"]
    for dist, res_id in results_dist:
        csv_lines.append(f"{dist:.2f},{res_id}")
    
    csv_path = output_dir / "residuos_raio.csv"
    csv_path.write_text("\n".join(csv_lines))
    
    # mfcc_residuos.txt: um resíduo por linha (sem distância)
    residues = [res_id for _, res_id in results_dist]
    mfcc_path = output_dir / "mfcc_residuos.txt"
    mfcc_path.write_text("\n".join(residues))
    
    return csv_path, mfcc_path
