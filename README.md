# Projet-Angle

## 📌 Description
Ce projet permet d'extraire et d'analyser les angles dièdres d'une structure moléculaire à partir d'un fichier **GROMACS (.gro)**. Il génère notamment un **diagramme de Ramachandran** pour visualiser les angles phi et psi.
Il inclut également l'extraction des frames d'une trajectoire **GROMACS (.xtc)** et la génération des diagrammes de Ramachandran pour chaque frame.  

## 📂 Structure du projet
``` 
projet-angle/
│── data/                 # contient les fichiers de données .gro, .pdb et .xtc
│ ├── md.gro
│ ├── start.gro
| ├── md_OK_dt100.xtc
│ ├── md.pdb
│ ├── start.pdb
│── scripts/              # contient les scripts python et shell
│ ├── pycache/
│ ├── main.py             # Script principal pour l'analyse
│ ├── packangle.py        # Bibliothèque des fonctions d'analyse
│ ├── run.sh              # Script bash pour exécuter les conversions avec GROMACS
│ ├── frame_extraction.py # Script pour extraire les frames en .gro à partir du fichier .xtc
│ ├── main_frames.py      # Script pour générer les diagrammes de Ramachandran des frames
│ ├── biopython_check.py  # Script pour analyser les angles phi et psi des fichiers PDB (BIO)
│── README.md             # Documentation du projet
```

####  Dossiers et fichiers :
- **data/** : Contient les fichiers `.gro` utilisés pour l'analyse.
- **scripts/** :
  - `main.py` : Programme principal pour l'extraction et l'affichage des angles dièdres.
  - `packangle.py` : Module contenant les fonctions de parsing et de calcul des angles.
  - `run.sh` : Script shell pour la conversion des fichiers de trajectoire GROMACS.
  - `frame_extraction.py` : Extraction des frames `.gro` à partir du fichier `.xtc`.
  - `main_frames.py` : Génération des diagrammes de Ramachandran pour chaque frame.
  - `biopython_check.py` : Analyse des fichiers `.pdb` avec Biopython pour extraire les angles phi et psi.
- **README.md** : Documentation du projet.

## 🚀 Installation et Pré-requis
Ce projet nécessite **Python 3** et les bibliothèques suivantes :
- `numpy`
- `matplotlib`
- `MDAnalysis`
- `Biopython`
-  `os`

Vous pouvez les installer avec :
```bash
pip install numpy matplotlib MDAnalysis biopython os
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
4- **Extraction des frames de la trajectoire .xtc** : 
Le script `frame_extraction.py` extrait 3001 frames en `.gro` à partir du fichier `.xtc`  :
```
python scripts/frame_extraction.py
```
cela génère les fichiers:
```
output/frames/frame_0000.gro
output/frames/frame_0001.gro
...
output/frames/frame_3000.gro
```
5- **Génération des diagrammes de Ramachandran pour chaque frame** : 
Le script `main_frames.py` calcule les angles phi et psi pour chaque frame et génère les diagrammes de Ramachandran correspondants : 
```
python scripts/main_frames.py
```
Les fichiers de sortie seront enregistrés dans `output/frames_ramachandran/`.

6- **Analyse des fichiers .pdb** : 
Le script `biopython_check.py` extrait et trace les diagrammes de Ramachandran pour les fichiers PDB :
```
python scripts/biopython_check.py
```
Les fichiers générés seront dans `biopython_ramachandran/`.
## 📬 Contacts
- Romain Gérard : romain.gerard@etu.u-paris.fr
- Gaelle Loutfi : gaelle.loutfi@etu.u-paris.fr
