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


#def compute_dihedral():

if __name__ == "__main__":
    print("Uniquement un package de fonctions utilisées dans main.py. Ne fais rien tout seul.")
