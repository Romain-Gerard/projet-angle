# TO BE UPDATED : Projet-Angle

## 📌 Description
Ce projet permet d'extraire et d'analyser les angles dièdres d'une structure moléculaire à partir d'un fichier **GROMACS (.gro)**. Il génère notamment un **diagramme de Ramachandran** pour visualiser les angles phi et psi.

## 📂 Structure du projet
``` 
projet-angle/
│── data/            # contient les fichiers de données .gro et .xtc
│ ├── md.gro
│ ├── start.gro
| ├── md_OK_dt100.xtc
│── scripts/         # contient les scripts python et shell
│ ├── pycache/
│ ├── main.py        # Script principal pour l'analyse
│ ├── packangle.py   # Bibliothèque des fonctions d'analyse
│ ├── run.sh         # Script bash pour exécuter les conversions avec GROMACS
│── README.md
```

####  Dossiers et fichiers :
- **data/** : Contient les fichiers `.gro` utilisés pour l'analyse.
- **scripts/** :
  - `main.py` : Programme principal pour l'extraction et l'affichage des angles dièdres.
  - `packangle.py` : Module contenant les fonctions de parsing et de calcul des angles.
  - `run.sh` : Script shell pour la conversion des fichiers de trajectoire GROMACS.
- **README.md** : Documentation du projet.

## 🚀 Installation et Pré-requis
Ce projet nécessite **Python 3** et les bibliothèques suivantes :
- `numpy`
- `matplotlib`

Vous pouvez les installer avec :
```bash
pip install numpy matplotlib
```
## 🛠️ Utilisation:
1- **Extraction des coordonnées des atomes** : Le script `packangle.py` contient une fonction `parse_gro()` qui extrait les coordonnées atomiques des résidus du fichier `.gro`.
```
from packangle import parse_gro

dict_aa = parse_gro("../data/start.gro")
print(dict_aa)
```
2- **Calcul des angles Phi et Psi** : Le script `main.py` extrait les coordonnées des atomes nécessaires et calcule les angles dièdres phi et psi. \n
Exécuter le script :
```
python scripts/main.py
```
3- **Génération du diagramme de Ramachandran** : 
Le script `main.py` génère un diagramme de Ramachandran basé sur les angles phi et psi extraits :
```
plt.scatter(phi, psi)
plt.xlabel("Phi")
plt.ylabel("Psi")
plt.title("Diagramme de Ramachandran")
plt.show()
```
4- **Execution du script run.sh** : 
Le script shell permet de convertir un fichier `.xtc` en `.pdb` en utilisant `gmx trjconv` :
```
bash scripts/run.sh
```
## 📬 Contacts
- Romain Gérard : romain.gerard@etu.u-paris.fr
- Gaelle Loutfi : gaelle.loutfi@etu.u-paris.fr
