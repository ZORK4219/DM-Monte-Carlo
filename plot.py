import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis le fichier
#data = np.loadtxt('Positions_xy.txt')

# Initialiser les listes pour les coordonnées
x_coords = []
y_coords = []

# Lire le fichier ligne par ligne
with open('Positions_xy_end.txt', 'r') as file:
    for line in file:
        # Si la ligne commence par 'A', c'est une ligne de données valide
        if line.startswith('A'):
            # Extraire les coordonnées x et y après le mot "A"
            parts = line.split()
            x_coords.append(float(parts[1]))  # Coordonnée x
            y_coords.append(float(parts[2]))  # Coordonnée y

# Tracer les coordonnées sur un graphique 2D
plt.figure(figsize=(8, 8))
plt.scatter(x_coords, y_coords, color='blue', marker='o', s=10)

# Ajouter des labels et un titre
plt.xlabel('Coordonnée x', fontsize=14)
plt.ylabel('Coordonnée y', fontsize=14)
#plt.title('Positions des partciles')

# Afficher le graphique
#plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Pour garder le même rapport entre les axes X et Y
plt.show()