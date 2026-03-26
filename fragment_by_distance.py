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
    
    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        
        record = line[0:6].strip()
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain_id = line[21].strip()
        res_seq = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        
        # Filtrar por cadeia
        if chain and chain_id != chain:
            continue
        
        # Ligante
        if record == "HETATM" and resname == ligand_resname:
            ligand_atoms.append((x, y, z))
        
        # Proteína
        elif record == "ATOM":
            res_id = f"{resname}{res_seq}"
            if res_id not in residues:
                residues[res_id] = []
            residues[res_id].append((x, y, z))
    
    return {"ligand_atoms": ligand_atoms, "residues": residues}


def centroid(coords: list) -> tuple:
    """Calcula centróide de uma lista de coordenadas."""
    n = len(coords)
    if n == 0:
        return (0, 0, 0)
    x = sum(c[0] for c in coords) / n
    y = sum(c[1] for c in coords) / n
    z = sum(c[2] for c in coords) / n
    return (x, y, z)


def distance_3d(p1: tuple, p2: tuple) -> float:
    """Distância euclidiana entre dois pontos 3D."""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)


def find_residues_in_cutoff(
    pdb_file: str,
    ligand_resname: str,
    chain: Optional[str] = None,
    cutoff: float = 8.0
) -> list:
    """
    Encontra resíduos dentro de cutoff Ångströms do ligante.
    
    Retorna:
        lista de tuplas (distância, resíduo_id), ordenada por distância
    """
    data = parse_pdb_ligand(pdb_file, ligand_resname, chain)
    
    if not data["ligand_atoms"]:
        return []
    
    ligand_centroid = centroid(data["ligand_atoms"])
    results = []
    
    for res_id, coords in data["residues"].items():
        res_centroid = centroid(coords)
        dist = distance_3d(ligand_centroid, res_centroid)
        
        if dist <= cutoff:
            results.append((dist, res_id))
    
    # Ordenar por distância
    results.sort(key=lambda x: x[0])
    
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
