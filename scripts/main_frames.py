import numpy as np
import matplotlib.pyplot as plt
import os
import packangle as pa

if __name__ == "__main__":
    # Définition du dossier de sortie pour stocker les diagrammes de Ramachandran
    output_folder = "../output/frames_ramachandran"
    os.makedirs(output_folder, exist_ok=True)  # Création du dossier s'il n'existe pas

    # Définition du dossier contenant les fichiers de frames .gro
    frames_folder = "../output/frames"

    # Boucle à travers toutes les frames de la trajectoire
    for frame_num in range(3001):  # On suppose que les frames sont numérotées de 0000 à 3000
        frame_filename = f"frame_{frame_num:04d}.gro"  # Format du nom du fichier
        frame_path = os.path.join(frames_folder, frame_filename)  # Construction du chemin d'accès au fichier
        
        # Vérification si le fichier existe
        if not os.path.exists(frame_path):
            print(f"Fichier non trouvé : {frame_path}")  # Affichage d'un message si le fichier est absent
            continue  # On passe à l'itération suivante
        
        print(f"Traitement du fichier : {frame_filename}")  # Indication de l'avancement

        # Chargement du fichier .gro et extraction des coordonnées atomiques
        dict_aa_frame = pa.parse_gro(frame_path)
        dict_phi_psi_coors_frame = pa.get_phi_psi_coors(dict_aa_frame)  # Extraction des coordonnées nécessaires
        dict_phi_psi_angles_frame, list_phi_psi_angles_frame = pa.get_phi_psi_angles(dict_phi_psi_coors_frame)  # Calcul des angles

        # Extraction des valeurs des angles phi et psi
        phi_frame = list_phi_psi_angles_frame[0]
        psi_frame = list_phi_psi_angles_frame[1]

        # Création du diagramme de Ramachandran
        plt.figure()
        plt.scatter(phi_frame, psi_frame, alpha=0.7, s=10)  # Tracé des points
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.xlabel("Phi (°)")
        plt.ylabel("Psi (°)")
        plt.title(f"Diagramme de Ramachandran - {frame_filename}")

        # Sauvegarde du diagramme dans le dossier spécifié
        output_path = os.path.join(output_folder, f"ramachandran_{frame_filename}.png")
        plt.savefig(output_path, dpi=300)  # Enregistrement avec une résolution élevée
        plt.close()  # Fermeture de la figure pour éviter les fuites de mémoire

        print(f"Diagramme sauvegardé : {output_path}")  # Confirmation de l'enregistrement
