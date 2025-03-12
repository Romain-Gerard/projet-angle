import numpy as np
import matplotlib.pyplot as plt
import os
import packangle as pa

if __name__ == "__main__":
    # Création du dossier pour stocker les diagrammes
    output_folder = "frames_ramachandran"
    os.makedirs(output_folder, exist_ok=True)

    # Génération des diagrammes de Ramachandran pour chaque frame
    for frame_num in range(3001):
        frame_filename = f"frame_{frame_num:04d}.gro"
        frame_path = frame_filename  # Les fichiers sont dans le même répertoire que le script
        
        # Vérification si le fichier existe
        if not os.path.exists(frame_path):
            print(f"Fichier non trouvé : {frame_path}")
            continue
        
        print(f"Traitement du fichier : {frame_filename}")
        dict_aa_frame = pa.parse_gro(frame_path)
        dict_phi_psi_coors_frame = pa.get_phi_psi_coors(dict_aa_frame)
        dict_phi_psi_angles_frame, list_phi_psi_angles_frame = pa.get_phi_psi_angles(dict_phi_psi_coors_frame)

        # Création du diagramme
        phi_frame = list_phi_psi_angles_frame[0]
        psi_frame = list_phi_psi_angles_frame[1]
        plt.figure()
        plt.scatter(phi_frame, psi_frame)
        plt.xlabel("Phi")
        plt.ylabel("Psi")
        plt.title(f"Diagramme de Ramachandran - {frame_filename}")
        
        # Sauvegarde du diagramme
        output_path = os.path.join(output_folder, f"ramachandran_{frame_filename}.png")
        plt.savefig(output_path)
        plt.close()

        print(f"Diagramme sauvegardé : {output_path}")
