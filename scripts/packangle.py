import numpy as np


def parse_gro(filename):
    """
    Lit et parse un fichier au format GROMACS (.gro) afin d’extraire les coordonnées
    atomiques de chaque résidu (à l’exception des solvants et des ions).

    Arguments:
        filename (str) :
            Le chemin vers le fichier .gro à analyser.

    Retourne:
        dict :
            Un dictionnaire imbriqué :
            
            {
                res_id_1: {
                    atom_name_1: (x1, y1, z1),
                    atom_name_2: (x2, y2, z2),
                    ...
                },
                res_id_2: {
                    atom_name_1: (x3, y3, z3),
                    atom_name_2: (x4, y4, z4),
                    ...
                },
                ...
            }

            • res_id_x (int) : Identifiant séquentiel unique du résidu.
            • atom_name_y (str) : Nom de l’atome.
            • (x, y, z) (tuple de float) : Coordonnées cartésiennes de l’atome.

    Description du fonctionnement:
        1. Le fichier est lu ligne par ligne (sauf les deux premières et la dernière
           qui ne contiennent pas d’informations utiles pour les résidus).
        2. Les lignes correspondant aux atomes de solvants (SOL) et d’ions (NA, CL)
           sont ignorées.
        3. Chaque nouvelle occurrence d’un nom de résidu (res_name) incrémente un
           identifiant numérique (res_id). Celui-ci sert de clé dans le dictionnaire
           principal.
        4. Les coordonnées (coor_x, coor_y, coor_z) de chaque atome sont extraites
           et enregistrées dans le dictionnaire associé au résidu en cours, avec
           pour clé le nom d’atome.

    Exemple d’utilisation:
        >>> dict_residus = parse_gro("ma_simulation.gro")
        >>> for res_id, atomes in dict_residus.items():
        ...     print(f"Résidu {res_id} : {atomes}")

    Remarque:
        - Les résidus de type SOL, NA et CL sont systématiquement ignorés.
        - L’identifiant de résidu (res_id) est purement séquentiel et ne correspond
          pas nécessairement à l’index de résidu réel dans le fichier .gro original.
    """
    dict_aa = {}
    with open(filename, "r") as filin:
        lines = filin.readlines()
        res_id = 0
        list_res_name = []

        # On ignore la première et la deuxième ligne (lines[0:2])
        # et la dernière ligne (lines[-1]), qui contient la taille de la boîte
        for line in lines[2:-1]:
            # On ne traite pas les atomes du solvant ou des ions
            if line[5:8].strip() in ["SOL", "NA", "CL"]:
                continue

            res_name = line[0:8].strip()
            atom_name = line[8:15].strip()
            coor_x = float(line[20:28])
            coor_y = float(line[28:36])
            coor_z = float(line[36:44])

            # Si on rencontre pour la première fois le résidu, on incrémente l'ID
            if res_name not in list_res_name:
                res_id += 1
                dict_aa[res_id] = {}  # chaque résidu a son propre dict
                list_res_name.append(res_name)

            # On stocke la position de chaque atome
            dict_aa[res_id][atom_name] = (coor_x, coor_y, coor_z)

    return dict_aa

def get_phi_psi_coors(dict_aa):
    """
    Pour chaque résidu entre le 2e et l'avant-dernier, extrait les coordonnées des atomes
    C précédent, N, CA, C et N suivant pour le calcul des angles phi et psi.

    Arguments:
        dict_aa (dict) :
            Un dictionnaire de coordonnées atomiques par résidu généré par la fonction parse_gro.
    
    Retourne:
        dict :
            Un dictionnaire de matrices numpy de coordonnées des atomes C précédent, N, CA, C et N suivant
            de chaque résidu. Les clés sont les mêmes que celles du dictionnaire d'entrée dict_aa.
            Représentation : {res_id: np.array([[x1, y1, z1],     # C précédent
                                                [x2, y2, z2],     # N
                                                [x3, y3, z3],     # CA
                                                [x4, y4, z4],     # C
                                                [x5, y5, z5]])},  # N suivant
                                                
                              res_id+1: np.array([[x6, y6, z6],   # C précédent
                                                  [x7, y7, z7],   # N
                                                  [x8, y8, z8],   # CA
                                                  [x9, y9, z9],   # C
                                                  [x10, y10, z10]])},  # N suivant
                              ...}
            Pour le premier et le dernier résidu, on n'a pas d'atome C précédent et N suivant respectivement.
        
    Fonctionnement:
        Pour chaque résidu, on vérifie si les atomes C, N et CA sont présents. Si c'est le cas,
        on stocke les coordonnées des atomes C précédent, N, CA, C et N suivant dans un tableau numpy.
        Pour le premier et le dernier résidu, on n'a pas d'atome C précédent et N suivant respectivement.

    Exemple d'utilisation:
        >>> dict_aa = parse_gro("ma_simulation.gro")
        >>> dict_phi_psi_coors = get_phi_psi_coors(dict_aa)
        >>> for res_id, mat_coors in dict_phi_psi_coors.items():
        ...     print(f"Résidu {res_id} : {mat_coors}")
    """
    dict_phi_psi = {}
    for res_id, atoms in dict_aa.items():
        if "C" in atoms and "N" in atoms and "CA" in atoms:
            # Si le résidu est le premier, on n'a pas d'atome C précédent,
            # seulement ceux pour le calcul de psi
            if res_id == 1:
                dict_phi_psi[res_id] = np.array([
                    atoms["N"],
                    atoms["CA"],
                    atoms["C"],
                    dict_aa[res_id + 1]["N"]
                ])
            # Si on est entre le 2e et l'avant-dernier résidu, on a tous les
            # atomes nécessaires pou rle calcul de phi et psi
            elif res_id > 1 and res_id < len(dict_aa):
                dict_phi_psi[res_id] = np.array([
                    dict_aa[res_id - 1]["C"],
                    atoms["N"],
                    atoms["CA"],
                    atoms["C"],
                    dict_aa[res_id + 1]["N"]
                ])
            # Si le résidu est le dernier, on n'a pas d'atome N suivant,
            # seulement ceux pour le calcul de phi
            elif res_id == len(dict_aa):
                dict_phi_psi[res_id] = np.array([
                    dict_aa[res_id - 1]["C"],
                    atoms["N"],
                    atoms["CA"],
                    atoms["C"]
                ])
    return dict_phi_psi

def get_khi1_khi2_coors(dict_aa):
    """
    Extrait les coordonnées des atomes N, CA, CB, CG, et CD de chaque résidu pour
    le calcul des angles khi1 et khi2.

    Arguments:
        dict_aa (dict) :
            Un dictionnaire de coordonnées atomiques par résidu généré par la fonction parse_gro.
        
    Retourne:
        dict :
            Un dictionnaire de matrices numpy de coordonnées des atomes N, CA, CB, CG et CD
            de chaque résidu. Les clés sont les mêmes que celles du dictionnaire d'entrée dict_aa.
    
    Fonctionnement:
        Pour chaque résidu, on vérifie si les atomes N, CA, CB, CG et CD sont présents.
        Si c'est le cas, on stocke les coordonnées des atomes dans un tableau numpy.

    Exemple d'utilisation:
        >>> dict_aa = parse_gro("ma_simulation.gro")
        >>> dict_khi_coors = get_khi1_khi2_coors(dict_aa)
        >>> for res_id, mat_coors in dict_khi_coors.items():
        ...     print(f"Résidu {res_id} : {mat_coors}")
    """
    dict_khi = {}
    for res_id, atoms in dict_aa.items():
        if "CB" in atoms and "CG" in atoms and "CD" in atoms:
            dict_khi[res_id] = np.array([
                atoms["N"],
                atoms["CA"],
                atoms["CB"],
                atoms["CG"],
                atoms["CD"]
            ])
    return dict_khi

def compute_dihedral(mat_coors):
    """
    Calcule l'angle dièdre en degrés entre 4 atomes consécutifs à partir de leurs coordonnées cartésiennes.
    
    Arguments:
        mat_coors (matrice numpy) : Une matrice de coordonnées cartésiennes de 4 atomes consécutifs.
                                   Chaque ligne correspond à un atome et chaque colonne à une coordonnée.
    
    Retourne:
        float : La valeur de l'angle dièdre en degrés compris entre -180° et +180°.
    
    Fonctionnement:
        1. On construit les vecteurs entre les atomes.
        2. On calcule les normales aux plans formés par les vecteurs.
        3. On calcule l'angle entre les deux normales en tenant compte du signe.
    
    Exemple d'utilisation:
        >>> mat_coors = np.array([
        ...     [0, 0, 0],
        ...     [1, 0, 0],
        ...     [1, 1, 0],
        ...     [1, 1, 1]
        ... ])
        >>> angle = compute_dihedral(mat_coors)
        >>> print(angle)
    """
    # On construit les vecteurs entre les atomes
    vec1 = mat_coors[1] - mat_coors[0]
    vec2 = mat_coors[2] - mat_coors[1]
    vec3 = mat_coors[3] - mat_coors[2]
    
    # On calcule les normales aux plans
    norm1 = np.cross(vec1, vec2)
    norm2 = np.cross(vec2, vec3)
    
    # Normalisation des vecteurs normaux
    norm1 = norm1 / np.linalg.norm(norm1)
    norm2 = norm2 / np.linalg.norm(norm2)
    
    # Calcul du produit scalaire pour déterminer l'angle
    cos_angle = np.dot(norm1, norm2)
    
    # Pour déterminer le signe de l'angle on utilise le produit mixte
    # avec le vecteur central comme référence
    vec2_norm = vec2 / np.linalg.norm(vec2)
    sign = np.sign(np.dot(np.cross(norm1, norm2), vec2_norm))
    
    # Calcul de l'angle en radians puis conversion en degrés
    angle_rad = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    if sign < 0:
        angle_rad = -angle_rad
    
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg

def get_phi_psi_angles(dict_phi_psi_coors):
    """
    Calcule les angles phi et psi pour chaque résidu à partir des coordonnées des atomes
    C précédent, N, CA, C et N suivant.

    Arguments:
        dict_phi_psi_coors (dict) :
            Un dictionnaire de matrices numpy de coordonnées des atomes C précédent, N, CA, C et N suivant
            de chaque résidu.

    Retourne:
        dict :
            Un dictionnaire des angles phi et psi de chaque résidu. Les clés sont les mêmes que celles
            du dictionnaire d'entrée dict_phi_psi_coors.
            Représentation : {res_id: (phi, psi),
                              res_id+1: (phi, psi),
                              ...}
        list :
            Une liste de listes contenant les angles phi et psi pour chaque résidu.
            Représentation : [[phi1, phi2, phi3, ...],
                              [psi1, psi2, psi3, ...]]
            
            Pour le premier et le dernier résidu, on n'a que l'angle psi et phi respectivement.
    
    Fonctionnement:
        Pour chaque résidu, on calcule les angles phi et psi à partir des coordonnées des atomes
        C précédent, N, CA, C et N suivant.
        Pour le premier et le dernier résidu, on ne calcule que l'angle psi et phi respectivement.
    
    Exemple d'utilisation:
        >>> dict_aa = parse_gro("ma_simulation.gro")
        >>> dict_phi_psi_coors = get_phi_psi_coors(dict_aa)
        >>> dict_phi_psi_angles = get_phi_psi_angles(dict_phi_psi_coors)
        >>> for res_id, angles in dict_phi_psi_angles.items():
        ...     print(f"Résidu {res_id} : {angles}")
    """
    dict_phi_psi_angles = {}
    list_phi_psi = [[], []]
    for res_id, mat_coors in dict_phi_psi_coors.items():
        # Si le résidu est le premier, on ne calcule que l'angle psi
        if res_id == 1:
            phi = None
            psi = compute_dihedral(mat_coors)
        # Si le résidu est le dernier, on ne calcule que l'angle phi
        elif res_id == len(dict_phi_psi_coors):
            phi = compute_dihedral(mat_coors)
            psi = None
        # Si on est entre le 2e et l'avant-dernier résidu, on calcule phi et psi
        else:
            phi = compute_dihedral(mat_coors[:4])
            psi = compute_dihedral(mat_coors[1:])
        # On ajoute les angles au dictionnaire
        dict_phi_psi_angles[res_id] = (phi, psi)
        # On ajoute les angles à la liste pour le diagramme de Ramachandran
        list_phi_psi[0].append(phi)
        list_phi_psi[1].append(psi)

    return dict_phi_psi_angles, list_phi_psi

def get_khi1_khi2_angles(dict_khi_coors):
    """
    Calcule les angles khi1 et khi2 pour chaque résidu à partir des coordonnées des atomes
    N, CA, CB, CG et CD.

    Arguments:
        dict_khi_coors (dict) :
            Un dictionnaire de matrices numpy de coordonnées des atomes N, CA, CB, CG et CD
            de chaque résidu.

    Retourne:
        dict :
            Un dictionnaire des angles khi1 et khi2 de chaque résidu. Les clés sont les mêmes que celles
            du dictionnaire d'entrée dict_khi_coors.
            Représentation : {res_id: (khi1, khi2),
                              res_id+1: (khi1, khi2),
                              ...}
        list :
            Une liste de listes contenant les angles khi1 et khi2 pour chaque résidu.
            Représentation : [[khi1_1, khi1_2, khi1_3, ...],
                              [khi2_1, khi2_2, khi2_3, ...]}
   
    Fonctionnement:
        Pour chaque résidu, on calcule les angles khi1 et khi2 à partir des coordonnées des atomes
        N, CA, CB, CG et CD.
    
    Exemple d'utilisation:
        >>> dict_aa = parse_gro("ma_simulation.gro")
        >>> dict_khi_coors = get_khi1_khi2_coors(dict_aa)
        >>> dict_khi_angles = get_khi1_khi2_angles(dict_khi_coors)
        >>> for res_id, angles in dict_khi_angles.items():
        ...     print(f"Résidu {res_id} : {angles}")
    """
    dict_khi_angles = {}
    list_khi = [[], []]
    for res_id, mat_coors in dict_khi_coors.items():
        khi1 = compute_dihedral(mat_coors[:4])
        khi2 = compute_dihedral(mat_coors[1:])
        dict_khi_angles[res_id] = (khi1, khi2)
        list_khi[0].append(khi1)
        list_khi[1].append(khi2)

    return dict_khi_angles, list_khi

if __name__ == "__main__":
    print("Package de fonctions utilisées dans main.py. Ne fais rien tout seul.")
