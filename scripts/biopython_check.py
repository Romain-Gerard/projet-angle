import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB

# Liste des fichiers PDB à analyser
pdb_files = ["../pdb/start.pdb", "../pdb/md.pdb"]

# Dossier pour stocker les diagrammes
output_folder = "biopython_ramachandran"
os.makedirs(output_folder, exist_ok=True)

# Analyseur PDB
parser = PDB.PDBParser(QUIET=True)

def extract_phi_psi_angles(pdb_filename):
    """
    Extrait les angles phi et psi d'un fichier PDB.
    """
    structure = parser.get_structure(pdb_filename, pdb_filename)
    phi_angles = []
    psi_angles = []

    # Analyse des modèles et chaînes
    for model in structure:
        for chain in model:
            poly = PDB.Polypeptide.Polypeptide(chain)
            for phi_psi in poly.get_phi_psi_list():
                if None not in phi_psi:  # Vérifie que les angles sont valides
                    phi_angles.append(np.degrees(phi_psi[0]))
                    psi_angles.append(np.degrees(phi_psi[1]))

    return phi_angles, psi_angles

def plot_ramachandran(pdb_filename, phi_angles, psi_angles):
    """
    Génère un diagramme de Ramachandran et l'enregistre.
    """
    plt.figure(figsize=(6,6))
    plt.scatter(phi_angles, psi_angles, alpha=0.7, s=10)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.xlabel("Phi (°)")
    plt.ylabel("Psi (°)")
    plt.title(f"Diagramme de Ramachandran - {os.path.basename(pdb_filename)}")

    # Sauvegarde du diagramme
    output_path = os.path.join(output_folder, f"ramachandran_{os.path.basename(pdb_filename)}.png")
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Diagramme sauvegardé : {output_path}")

# Exécuter l'analyse sur chaque fichier PDB
for pdb_file in pdb_files:
    if os.path.exists(pdb_file):
        print(f"Traitement de {pdb_file}...")
        phi, psi = extract_phi_psi_angles(pdb_file)
        plot_ramachandran(pdb_file, phi, psi)
    else:
        print(f"Fichier non trouvé : {pdb_file}")
