import MDAnalysis as mda

# Charger la trajectoire en utilisant start.gro comme fichier de topologie
u = mda.Universe("../data/start.gro", "../data/md_OK_dt100.xtc")

# Parcourir toutes les frames de la trajectoire
for ts in u.trajectory:
    frame_number = ts.frame  # Récupérer l'index de la frame actuelle
    output_filename = f"frame_{frame_number:04d}.gro"  # Générer un nom de fichier, ex : frame_0001.gro

    # Écrire la frame actuelle dans un fichier .gro
    with mda.Writer(output_filename, n_atoms=u.atoms.n_atoms) as w:
        w.write(u.atoms)  # Sauvegarder les atomes de la frame actuelle

# Afficher le dernier fichier sauvegardé
print(f"Sauvegardé : {output_filename}")
