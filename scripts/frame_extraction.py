import MDAnalysis as mda
import os

# Charger la trajectoire en utilisant start.gro comme fichier de topologie
u = mda.Universe("../data/start.gro", "../data/md_OK_dt100.xtc")

# Définition du dossier de sortie
output_folder = "../output/frames"
os.makedirs(output_folder, exist_ok=True)  # Créer le dossier s'il n'existe pas

# Parcourir toutes les frame de la trajectoire
for ts in u.trajectory:
    frame_number = ts.frame  # Récupérer l'index de la frame actuelle
    output_filename = os.path.join(output_folder, f"frame_{frame_number:04d}.gro")  # Chemin complet du fichier

    # Écrire la frame actuelle dans un fichier .gro
    with mda.Writer(output_filename, n_atoms=u.atoms.n_atoms) as w:
        w.write(u.atoms)  # Sauvegarder les atomes de la frame actuelle

print(f"Frames sauvegardées dans : {output_folder}")
