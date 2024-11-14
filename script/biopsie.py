import numpy as np
import matplotlib.pyplot as plt

# CODE 1 - Extension à 15%
# Variables
ex = 0.08  # Déformation dans la direction x
ey = 0.024  # Déformation dans la direction y
exy = 0.024  # Déformation de cisaillement

# 1. Calcul du coefficient de Poisson
v = -(ey / ex)
# Le coefficient de Poisson est négatif, donc le matériau étudié est auxétique.

# 2. Calcul des directions principales et déformations principales
# Déformations principales
eX = (ex + ey) / 2 + 0.5 * np.sqrt((ex - ey) ** 2 + 4 * exy ** 2)
eY = (ex + ey) / 2 - 0.5 * np.sqrt((ex - ey) ** 2 + 4 * exy ** 2)

# Directions principales
PX = (180 / np.pi) * 0.5 * np.arctan(2 * exy / (ex - ey))
PY = PX + 90

# Calcul des vecteurs pour les déformations en fonction de phi
Phi = np.linspace(0, 360, 180)  # Angles de 0° à 360°
PhiRad = np.deg2rad(Phi)  # Conversion en radians

# Déformation normale et de cisaillement
en = (ex + ey) / 2 + ((ex - ey) / 2) * np.cos(2 * PhiRad) + exy * np.sin(2 * PhiRad)
yn = (ex - ey) / 2 * np.sin(2 * PhiRad) + exy * np.cos(2 * PhiRad)

# Tracé des déformations
plt.figure()
plt.polar(PhiRad, en, 'b', linewidth=2, label='Déformation ε')  # Déformation normale en bleu
plt.polar(PhiRad, yn, 'r--', linewidth=2, label='Déformation γ')  # Déformation de cisaillement en rouge pointillé
plt.title('Extension à 15%')
plt.legend()
plt.show()

# CODE 2 - Extension à 10%
# Variables
ex = 0.146  # Déformation dans la direction x
ey = 0.042  # Déformation dans la direction y
exy = 1.48  # Déformation de cisaillement

# 1. Calcul du coefficient de Poisson
v = -(ey / ex)
# Le coefficient de Poisson est négatif, donc le matériau étudié est auxétique.

# 2. Calcul des directions principales et déformations principales
# Déformations principales
eX = (ex + ey) / 2 + 0.5 * np.sqrt((ex - ey) ** 2 + 4 * exy ** 2)
eY = (ex + ey) / 2 - 0.5 * np.sqrt((ex - ey) ** 2 + 4 * exy ** 2)

# Directions principales
PX = (180 / np.pi) * 0.5 * np.arctan(2 * exy / (ex - ey))
PY = PX + 90

# Calcul des vecteurs pour les déformations en fonction de phi
Phi = np.linspace(0, 360, 180)  # Angles de 0° à 360°
PhiRad = np.deg2rad(Phi)  # Conversion en radians

# Déformation normale et de cisaillement
en = (ex + ey) / 2 + ((ex - ey) / 2) * np.cos(2 * PhiRad) + exy * np.sin(2 * PhiRad)
yn = (ex - ey) / 2 * np.sin(2 * PhiRad) + exy * np.cos(2 * PhiRad)

# Tracé des déformations
plt.figure()
plt.polar(PhiRad, en, 'b', linewidth=2, label='Déformation ε')  # Déformation normale en bleu
plt.polar(PhiRad, yn, 'r--', linewidth=2, label='Déformation γ')  # Déformation de cisaillement en rouge pointillé
plt.title('Extension à 10%')
plt.legend()
plt.show()
