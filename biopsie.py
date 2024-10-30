import matplotlib.pyplot as plt
import math
import numpy as np

# Deformations
ex10 = 0.146
ey10 = 0.042
exy10 = 1.48

ex15 = 0.08
ey15 = 0.024
exy15 = 0.024

# Valeurs de phi1 de 0 à 180 degrés avec un pas de 10
phi1_values = np.arange(0, 361, 10)
en10_values = []
en15_values = []
gamma10_values = []
gamma15_values = []

# Calcul de en10, en15, gamma10 et gamma15 pour chaque valeur de phi1
for phi1 in phi1_values:
    # Calcul de en10
    en10 = (ex10 + ey10) / 2 + ((ex10 - ey10) / 2) * math.cos(2 * math.radians(phi1)) + exy10 * math.sin(2 * math.radians(phi1))
    en10_values.append(en10)
    
    # Calcul de en15
    en15 = (ex15 + ey15) / 2 - ((ex15 - ey15) / 2) * math.cos(2 * math.radians(phi1)) - exy15 * math.sin(2 * math.radians(phi1))
    en15_values.append(en15)
    
    # Calcul de gamma10
    gamma10 = math.sin(2 * math.radians(phi1)) * (ex10 + ey10) + 2 * exy10 * math.cos(2 * math.radians(phi1))
    gamma10_values.append(gamma10)
    
    # Calcul de gamma15
    gamma15 = math.sin(2 * math.radians(phi1)) * (ex15 + ey15) + 2 * exy15 * math.cos(2 * math.radians(phi1))
    gamma15_values.append(gamma15)

# Création d'un graphique polaire
theta = np.radians(phi1_values)  # Convertir en radians pour le graphique polaire

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

# Tracer en10
ax.plot(theta, en10_values,label='en10', color='blue', linestyle='-')
# Tracer en15
ax.plot(theta, en15_values, label='en15', color='orange', linestyle='-')
# Tracer gamma10
ax.plot(theta, gamma10_values, label='gamma10', color='green', linestyle='--')
# Tracer gamma15
ax.plot(theta, gamma15_values, label='gamma15', color='red', linestyle='--')

# Ajouter des éléments au graphique
ax.set_title("Variation de en10, en15, gamma10 et gamma15 en fonction de phi1", va='bottom')
ax.set_xlabel("phi1 (degrés)")
ax.set_ylabel("Valeurs")
ax.legend()

plt.show()
