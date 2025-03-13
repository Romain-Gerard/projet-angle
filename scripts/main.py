import numpy as np
import matplotlib.pyplot as plt
import packangle as pa

if __name__ == "__main__":
    # On parse le fichier .gro
    dict_aa = pa.parse_gro("../data/start.gro")

    # On extrait les coordonnées des atomes nécessaires pour le calcul des angles phi et psi
    # Pour chaque résidu entre le 2e et l'avant-dernier, on a les coordonnées des
    # atomes C précédent, N, CA, C et N suivant.
    dict_phi_psi_coors= pa.get_phi_psi_coors(dict_aa)
    print(len(dict_phi_psi_coors))

    # On extrait les coordonnées des atomes nécessaires pour le calcul des angles khi1 et khi2
    # Pour chaque résidu, on a les coordonnées des atomes N, CA, CB, CG et CD s'ils existent.
    dict_khi1_khi2_coors = pa.get_khi1_khi2_coors(dict_aa)
    print(len(dict_khi1_khi2_coors))

    # On calcule les angles phi et psi
    dict_phi_psi_angles, list_phi_psi_angles = pa.get_phi_psi_angles(dict_phi_psi_coors)
    print(len(dict_phi_psi_angles))
    print(list_phi_psi_angles)
    
    # On affiche le diagramme de Ramachandran pour les angles phi et psi
    phi = list_phi_psi_angles[0]
    psi = list_phi_psi_angles[1]
    plt.figure(figsize=(8, 8))
    plt.scatter(phi, psi, alpha=0.7, edgecolor='k', s=50)
    plt.xlabel("φ°", fontsize=12)
    plt.ylabel("ψ°", fontsize=12)
    plt.title("Diagramme de Ramachandran", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.xticks(np.arange(-180, 181, 60))
    plt.yticks(np.arange(-180, 181, 60))
    plt.show
    #plt.savefig("diagramme_test.png")


