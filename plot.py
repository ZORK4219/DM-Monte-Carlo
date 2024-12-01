import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis le fichier
data = np.loadtxt('Positions_xy.txt')

# Extraire les coordonnées x et y
x = data[:, 0]
y = data[:, 1]

# Taille de la boîte
L = 45

# Créer une figure et un sous-graphique
fig, ax = plt.subplots()

# Tracer les particules
ax.scatter(x, y, c='r', marker='o', label='Particules')

# Dessiner la boîte
box_x = -L / 2
box_y = -L / 2
ax.add_patch(plt.Rectangle((box_x, box_y), L, L, fill=False, edgecolor='blue', linewidth=2, label='Boîte'))

# Ajouter des labels et une légende
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Réseau de particules avec boîte limite')
ax.legend()

# Ajuster les limites du graphique pour inclure toute la boîte
ax.set_xlim(box_x - 5, box_x + L + 5)
ax.set_ylim(box_y - 5, box_y + L + 5)

# Afficher le plot
plt.axis('equal')  # Assurer une échelle égale pour x et y
plt.show()
