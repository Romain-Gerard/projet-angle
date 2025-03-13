import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB

# Liste des fichiers PDB à analyser
pdb_files = ["../data/pdb/start.pdb", "../data/pdb/md.pdb"]

# Définition du dossier de sortie pour stocker les diagrammes
output_parent_folder = "../output"  # Dossier parent
output_folder = os.path.join(output_parent_folder, "biopython_ramachandran")  # Sous-dossier spécifique
os.makedirs(output_folder, exist_ok=True)  # Création du dossier si inexistant

# Initialisation de l'analyseur PDB de Biopython
parser = PDB.PDBParser(QUIET=True)

def extract_phi_psi_angles(pdb_filename):
    """
    Extrait les angles phi et psi d'une structure protéique contenue dans un fichier PDB.
    
    Arguments :
        pdb_filename (str) : Chemin du fichier PDB à analyser.
    
    Retourne :
        tuple : Deux listes contenant les valeurs des angles phi et psi (en degrés).
    """
    structure = parser.get_structure(pdb_filename, pdb_filename)  # Chargement de la structure
    phi_angles = []
    psi_angles = []

    # Analyse des modèles et chaînes de la structure
    for model in structure:
        for chain in model:
            poly = PDB.Polypeptide.Polypeptide(chain)  # Création d'un objet polypeptide pour extraction des angles
            for phi_psi in poly.get_phi_psi_list():
                if None not in phi_psi:  # Vérification que les angles ne sont pas vides
                    phi_angles.append(np.degrees(phi_psi[0]))  # Conversion en degrés et ajout à la liste
                    psi_angles.append(np.degrees(phi_psi[1]))
    
    return phi_angles, psi_angles

def plot_ramachandran(pdb_filename, phi_angles, psi_angles):
    """
    Génère un diagramme de Ramachandran et enregistre l'image correspondante.
    
    Arguments :
        pdb_filename (str) : Nom du fichier PDB analysé.
        phi_angles (list) : Liste des angles phi extraits.
        psi_angles (list) : Liste des angles psi extraits.
    """
    plt.figure(figsize=(6,6))
    plt.scatter(phi_angles, psi_angles, alpha=0.7, s=10)  # Tracé des angles phi/psi sous forme de points
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.xticks(np.arange(-180, 181, 60))  # Ajouter des graduations sur l'axe X
    plt.yticks(np.arange(-180, 181, 60))  # Ajouter des graduations sur l'axe Y
    plt.xlabel("Phi (°)")
    plt.ylabel("Psi (°)")
    plt.title(f"Diagramme de Ramachandran - {os.path.basename(pdb_filename)}")
    
    # Définition du chemin de sauvegarde
    output_path = os.path.join(output_folder, f"ramachandran_{os.path.basename(pdb_filename)}.png")
    plt.savefig(output_path, dpi=300)  # Sauvegarde de l'image en haute résolution
    plt.close()  # Fermeture de la figure pour éviter une surcharge mémoire

    print(f"Diagramme sauvegardé : {output_path}")

# Parcours et analyse des fichiers PDB fournis
for pdb_file in pdb_files:
    if os.path.exists(pdb_file):  # Vérification de l'existence du fichier
        print(f"Traitement du fichier : {pdb_file}...")
        phi, psi = extract_phi_psi_angles(pdb_file)  # Extraction des angles
        plot_ramachandran(pdb_file, phi, psi)  # Génération et sauvegarde du diagramme
    else:
        print(f"Fichier non trouvé : {pdb_file}")

