#!/usr/bin/env python3
"""
generate_database.py
====================
Parser de PDB para gerar DATABASE.txt com informações de átomos,
cargas e tipos atômicos (para uso no pipeline MFCC).
"""

from pathlib import Path
from typing import Optional
import re


def format_resname_4digits(resname: str, resid: int) -> str:
    """Formata RES + número com exatamente 4 dígitos: ALA0027, GLU0104."""
    return f"{resname.strip()}{resid:04d}"


def parse_pdb(pdb_file: str, chain: Optional[str] = None, ligand_resname: str = "LIG") -> tuple:
    """
    Parse um arquivo PDB e extrai informações de átomos.
    
    Retorna:
        (db_lines, res_list)
        - db_lines: lista de tuplas (residue_id ou None para ligante, atom_info_dict)
        - res_list: lista ordenada de resíduos proteicos únicos
    """
    pdb_lines = Path(pdb_file).read_text().splitlines()
    
    db_lines = []
    res_set = set()
    
    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        
        record   = line[0:6].strip()           # ATOM or HETATM
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        resname  = line[17:20].strip()
        chain_id = line[21].strip()
        res_seq  = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        elem = line[76:78].strip() or atom_name[0]  # Fallback: primeiro caractere do nome
        
        # Filtrar por cadeia
        if chain and chain_id != chain:
            continue
        
        # Separar ligante de proteína
        if record == "HETATM" and resname == ligand_resname:
            res_id = None  # Ligante
        else:
            res_id = format_resname_4digits(resname, res_seq)
            if record == "ATOM":
                res_set.add(res_id)
        
        # Estimar carga (simplificado)
        charge = estimate_charge(atom_name, resname)
        
        atom_info = {
            "atom_num": atom_num,
            "atom_name": atom_name,
            "resname": resname,
            "res_seq": res_seq,
            "x": x,
            "y": y,
            "z": z,
            "elem": elem,
            "charge": charge,
        }
        
        db_lines.append((res_id, atom_info))
    
    res_list = sorted(list(res_set))
    return db_lines, res_list


def estimate_charge(atom_name: str, resname: str) -> float:
    """
    Estimativa simplificada de carga parcial.
    Valores reais viriam de AMBER/GAFF/etc.
    """
    # Átomos negativos
    if atom_name in {"OD1", "OD2", "OE1", "OE2"}:  # ASP, GLU
        return -0.5
    if atom_name == "O" and resname in {"ASP", "GLU"}:
        return -0.5
    
    # Átomos positivos
    if atom_name in {"NZ", "NH1", "NH2"}:  # LYS, ARG
        return 0.5
    if atom_name == "N" and resname in {"LYS", "ARG"}:
        return 0.3
    
    # Átomos polares neutros
    if atom_name in {"O", "OG", "OH"}:
        return -0.3
    if atom_name in {"N", "ND", "NE"}:
        return 0.1
    
    return 0.0


def write_outputs(db_lines: list, res_list: list, output_dir: Path) -> tuple:
    """
    Escreve DATABASE.txt e lista_residuo.txt.
    
    Retorna:
        (database_path, residue_list_path)
    """
    output_dir = Path(output_dir)
    
    # DATABASE.txt: cada linha = RESNAME TAB ELEM x y z TAB CHARGE TAB ATOM_TYPE
    db_lines_text = []
    for res_id, atom_info in db_lines:
        display_res = atom_info["resname"]
        elem = atom_info["elem"]
        x, y, z = atom_info["x"], atom_info["y"], atom_info["z"]
        charge = atom_info["charge"]
        atom_type = atom_info["atom_name"]
        
        line = f"{display_res}\t{elem} {x:.4f} {y:.4f} {z:.4f}\t{charge:.4f}\t{atom_type}"
        db_lines_text.append(line)
    
    db_path = output_dir / "DATABASE.txt"
    db_path.write_text("\n".join(db_lines_text))
    
    # lista_residuo.txt: um resíduo por linha
    lr_path = output_dir / "lista_residuo.txt"
    lr_path.write_text("\n".join(res_list))
    
    return db_path, lr_path
