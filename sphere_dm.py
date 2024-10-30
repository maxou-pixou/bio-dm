import numpy as np
import matplotlib.pyplot as plt
 
# Données des matériaux
E_acier = 210e9  # Module de Young de l'acier (Pa)
nu_acier = 0.3   # Coefficient de Poisson de l'acier
E_poly = 10e6    # Module de Young du polymère (Pa)
nu_poly = 0.45   # Coefficient de Poisson du polymère
 
# Paramètres pour l'analyse
rayon_bille = 0.01  # Rayon de la bille (10 mm)
force = 50          # Force appliquée (50 N)
 
# Calcul du module d'Young équivalent pour acier-polymère
def calcul_E_equivalent(E1, nu1, E2, nu2):
    return 1 / ((1 - nu1**2) / E1 + (1 - nu2**2) / E2)
 
# Fonction pour calculer la pression maximale p0
def pression_maximale(F, R, E_eq):
    a = ((3 * F * R) / (4 * E_eq))**(1/3)  # Rayon de contact
    p0 = (3 * F) / (2 * np.pi * a**2)       # Pression maximale
    return p0, a
 
# Fonctions pour calculer les contraintes
def contraintes_interieures_r(r, a, nu, p0):
    term1 = (1 - 2 * nu) / 3
    term2 = (a**2 / r**2)
    return p0 * (term1 * term2 * (1 - (1 - (r**2 / a**2))**(3/2)) - np.sqrt(1 - (r**2 / a**2)))
 
def contraintes_interieures_theta(r, a, nu, p0):
    term1 = -(1 - 2 * nu) / 3
    term2 = (a**2 / r**2)
    return p0 * (term1 * term2 * (1 - (1 - (r**2 / a**2))**(3/2)) - 2 * nu * np.sqrt(1 - (r**2 / a**2)))
 
def contraintes_axiales(r, a, p0):
    return -p0 * np.sqrt(1 - (r / a)**2)
 
def contraintes_exterieures_r(r, a, nu, p0):
    return p0 * (1 - 2 * nu) * (a**2) / (3 * r**2)
 
# Calcul du module de Young équivalent
E_eq_acier_poly = calcul_E_equivalent(E_acier, nu_acier, E_poly, nu_poly)
 
# Calcul de la pression maximale
p0, max_a = pression_maximale(force, rayon_bille, E_eq_acier_poly)
 
# Créer des valeurs de r de -1.5a à 1.5a
r_values = np.linspace(-1.5*max_a, 1.5*max_a, 1000)  # Rayon d'affichage
 
# Calcul des contraintes
sigma_r = np.zeros(len(r_values))
sigma_theta = np.zeros(len(r_values))
sigma_z = np.zeros(len(r_values))
 
for i, r in enumerate(r_values):
    if abs(r) <= max_a:
        sigma_r[i] = contraintes_interieures_r(abs(r), max_a, nu_poly, p0)
        sigma_theta[i] = contraintes_interieures_theta(abs(r), max_a, nu_poly, p0)
        sigma_z[i] = contraintes_axiales(abs(r), max_a, p0)
    else:
        sigma_r[i] = contraintes_exterieures_r(abs(r), max_a, nu_poly, p0)
 
# Création du graphique
plt.figure(figsize=(10, 6))
 
# Tracer les contraintes par rapport à r/a
plt.plot(r_values / max_a, sigma_r / p0, linestyle='-', label=r'$\sigma_r$', color='green')
plt.plot(r_values / max_a, sigma_theta / p0, linestyle='-', label=r'$\sigma_\theta$', color='blue')
plt.plot(r_values / max_a, sigma_z / p0, linestyle='-', label=r'$\sigma_z$', color='red')
 
# Personnalisation des axes et de l'apparence
plt.axhline(0, color='black', lw=1)  # Ligne horizontale au niveau 0
plt.axvline(0, color='black', lw=1)  # Ligne verticale au niveau 0
plt.axvline(1, color='gray', linestyle='--', lw=0.5)
plt.axvline(-1, color='gray', linestyle='--', lw=0.5)
 
# Labels et légende
plt.xlabel(r"$r/a$")
plt.ylabel(r"Contraintes normalisées $\sigma/p_0$")
plt.legend(loc="lower right")
plt.title("Contraintes à l'intérieur et à l'extérieur de la surface de contact")
 
# Inversion de l'axe des y pour suivre le modèle de référence
# plt.gca().invert_yaxis()
 
# Ajustement des marges pour éviter les chevauchements
plt.tight_layout()
plt.grid(True)
 
# Affichage
plt.show()