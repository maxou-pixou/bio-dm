import numpy as np
import matplotlib.pyplot as plt

# Paramètres physiques
F_values = [50, 100]  # Forces appliquées (en Newtons)
radii = [10e-3, 5e-3, 2e-3, 500e-6]  # Différents rayons de la sphère (en mètres)
E_sphere = 210e9  # Module de Young de la sphère (acier, en Pascals)
nu_sphere = 0.3  # Coefficient de Poisson de la sphère
E_polymer_high = 210e9  # Module de Young du polymère (en Pascals)
nu_polymer = 0.45  # Coefficient de Poisson du polymère
E_polymer_low = 10e6  # Nouveau module de Young du polymère (en Pascals)

# Fonction pour calculer E* (module d'élasticité réduit)
def reduced_modulus(E1, nu1, E2, nu2):
    return 1 / ((1 - nu1**2) / E1 + (1 - nu2**2) / E2)

# Fonction pour calculer le rayon de la zone de contact
def contact_radius(F, R, E_star):
    return ((3 * F * R) / (4 * E_star))**(1/3)

# Fonction pour calculer la pression de contact maximale
def max_pressure(F, a):
    return (3 * F) / (2 * np.pi * a**2)

# Fonction pour calculer les contraintes radiales, circonférentielles et axiales
def contraintes_polaires(r, a, p0, nu):
    # Contraintes radiales σ_r
    sigma_r = p0 * (1 - 2 * nu / 3) * ((a / r)**2) * \
              (1 - (1 - (r / a)**2)**(3/2)) - p0 * (1 - (r / a)**2)**(1/2)
    
    # Contraintes circonférentielles σ_θ
    sigma_theta = -p0 * (1 - 2 * nu / 3) * ((a / r)**2) * \
                  (1 - (1 - (r / a)**2)**(3/2)) - 2 * nu * p0 * (1 - (r / a)**2)**(1/2)
    
    # Contraintes axiales σ_z
    sigma_z = -p0 * (1 - (r / a)**2)**(1/2)
    
    return sigma_r, sigma_theta, sigma_z

# Calcul du module d'élasticité réduit pour les deux valeurs de E_polymer
E_star_high = reduced_modulus(E_sphere, nu_sphere, E_polymer_high, nu_polymer)
E_star_low = reduced_modulus(E_sphere, nu_sphere, E_polymer_low, nu_polymer)

# Création des sous-graphiques pour les contraintes
fig, axs = plt.subplots(2, 2, figsize=(14, 14))

# Configuration des titres des sous-graphiques
titles = [
    'Contraintes radiales (σ_r)',
    'Contraintes circonférentielles (σ_θ)',
    'Contraintes axiales (σ_z)',
    'Comparaison σ_r, σ_θ, σ_z pour E = 10 MPa, F = 100 N'
]

# Boucle sur les forces et les modules d'élasticité
for i, F in enumerate(F_values):
    for j, (E_star, E_polymer_label) in enumerate([(E_star_high, '210 GPa'), (E_star_low, '10 MPa')]):
        for R in radii:
            # Calcul du rayon de la zone de contact
            a = contact_radius(F, R, E_star)
            
            # Rayon de la zone de contact
            r = np.linspace(0, a, 100)
            
            # Calcul de la pression de contact maximale
            p0 = max_pressure(F, a)
            
            # Calcul des contraintes polaires
            sigma_r, sigma_theta, sigma_z = contraintes_polaires(r, a, p0, nu_polymer)
            
            # Tracé des courbes polaires
            label = f'R = {R*1e3:.1f} mm, a = {a*1e3:.2f} mm'
            axs[0, 0].plot(r*1e3, sigma_r / 1e6, label=label)  # σ_r
            axs[0, 1].plot(r*1e3, sigma_theta / 1e6, label=label)  # σ_θ
            axs[1, 0].plot(r*1e3, sigma_z / 1e6, label=label)  # σ_z
            
# Configuration des sous-graphiques
for ax, title in zip(axs.flat, titles):
    ax.set_title(title)
    ax.set_xlabel('Rayon de la zone de contact (mm)')
    ax.set_ylabel('Contraintes (MPa)')
    ax.legend()
    ax.grid(True)

# Ajustement de l'espacement entre les sous-graphiques
plt.tight_layout(pad=5.0)

# Affichage du graphique
plt.show()
