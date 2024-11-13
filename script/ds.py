import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from mpl_toolkits.mplot3d import Axes3D
from sklearn.linear_model import LinearRegression

# Paramètres physiques
F_values = [50e-3]  # Forces appliquées (en Newtons)
radii= [2.5e-3] # Différents rayons de la sphère (en mètres)

E_sphere = 210e9  # Module de Young de la sphère (acier, en Pascals)
nu_sphere = 0.3  # Coefficient de Poisson de la sphère
E_polymer_high = 10e6  # Module de Young du polymère (en Pascals)
nu_polymer = 0.45  # Coefficient de Poisson du polymère
E_values =  [10e6, 210e9]


# Calcul du module d'Young équivalent pour acier-polymère
def reduced_modulus(E2, nu2):
    return  (1 - nu2**2) / E2
 
# Fonction pour calculer la pression maximale p0
def p_a(F, R,Z,E):
    a = np.sqrt(R*Z)
    p0 = ((6 * F * E**2) / ((np.pi)**3 * R**2))**1/3       # Pression maximale
    
    return p0,a



def r_equivalent(R1, R2):
    return np.reciprocal((1/R1) + (1/R2))


# Fonctions pour calculer les contraintes
def contr_int_r(r, a, p0, nu):
    term1 = (1 - 2 * nu) / 3
    term2 = (a**2 / r**2)
    return p0 * (term1 * term2 * (1 - (1 - (r**2 / a**2))**(3/2)) - np.sqrt(1 - (r**2 / a**2)))
 
def contr_int_theta(r, a, p0, nu):
    term1 = -(1 - 2 * nu) / 3
    term2 = (a**2 / r**2)
    return p0 * (term1 * term2 * (1 - (1 - (r**2 / a**2))**(3/2)) - 2 * nu * np.sqrt(1 - (r**2 / a**2)))
 
def contr_axi(r, a, p0):
    return -p0 * np.sqrt(1 - (r / a)**2)

def contr_ext_r(r, a, p0, nu):
    return p0 * (1 - 2 * nu) * (a**2) / (3 * r**2)

def contr_ext_theta(r, a, p0, nu):
    return - p0 * (1 - 2 * nu) * (a**2) / (3 * r**2)
 
def contr_cis_rt(z, a, p0, nu):
    return p0* (-(1+nu)*(1-(z/a)*np.arctan(a/z))+0.5*(1+(z/a)**2)**(-1))

def contr_cis_z(z, a):
    return (-(1+(z/a)**2)**(-1))
    
#################################################################################################################


    E_star = reduced_modulus(E_sphere, nu_sphere, E_sphere, nu_sphere)
    fig, ax = plt.subplots(figsize=(8, 6))
    num_radii = len(radii)
    colors = plt.cm.Accent(np.linspace(0, 1, num_radii))  # Utilisation de la colormap viridis

    # Initialiser le dictionnaire pour stocker les contraintes en fonction des rayons
    sigma_z_l = {}

    # Parcours de chaque rayon dans la liste radii pour calculer les contraintes
    for i, R in enumerate(radii):
        p0, a = p_a(F, R/2, E_star)
        
        # Générer les valeurs de r pour ce rayon spécifique
        r_values = np.linspace(-1.5 * a, 1.5 * a, 1000)
        
        # Calcul des contraintes pour chaque r
        sigma_z = np.zeros(len(r_values))
        for j, r in enumerate(r_values):
            if abs(r) <= a:
                sigma_z[j] = contr_axi(abs(r), a, p0)
            else : 
                sigma_z[j] = None
        
        # Ajouter les contraintes calculées pour chaque R au dictionnaire
        sigma_z_l[i] = (r_values, sigma_z)  # Associe les valeurs de r et les contraintes pour chaque rayon

    # Boucle pour afficher chaque courbe stockée dans sigma_r_l
    for i, (r_values, sigma_z) in sigma_z_l.items():
        R = radii[i]
        a = np.max(np.abs(r_values)) / 1.5  # Calcule a à partir de r_values pour chaque rayon R
        # Tracé des contraintes et des contraintes symétrisées avec la même couleur
        color = colors[i]
        # Tracer les contraintes pour chaque rayon avec une étiquette unique
        label = f'R = {R*1e3/2:.1f}* mm, R = {R*1e3:.1f} mm, F = {F} N, Sphère-Sphère'  # Label en mm pour R
        ax.plot(r_values / a, sigma_z / p0, label=f'σ_z ({label})', color=color)
        ax.plot(-r_values / a, sigma_z / p0, label='_nolegend_', color=color)  # Courbe symétrique
      

    # Configuration des labels et du titre
    ax.set_title('Contraintes axiales en fonction du rayon pour une sphère-sphère et F = 50N')
    ax.set_xlabel('Rayon de la zone de contact (mm)')
    ax.set_ylabel('Contraintes (MPa)')
    ax.legend()
    ax.grid(True)
    ax.invert_yaxis()
    ax.set_xlim(-1.5, 1.5)

    # Afficher le graphique dans le frame Tkinter
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
def plot_rayon_z_sppl(F, frame):

    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(8, 6))
    num_radii = len(radii)
    

    # Initialiser le dictionnaire pour stocker les contraintes en fonction des rayons
    sigma_z_l = {}
    max_d = 0.000637
    max_d1 = 0.001743
    # Parcours de chaque rayon dans la liste radii pour calculer les contraintes
    for i, R in enumerate(radii):
        # Générer les valeurs de r pour ce rayon spécifique
        r_values = np.linspace(0,max_d, 1000)
        r_values1 = np.linspace(0,max_d1, 1000)

        p0, a = p_a(50e-3, R, r_values, 43000)
        p01,a1 = p_a(50e-3, R, r_values1, 9000)
        
        sigma_z=contr_axi(r_values, a, p0)
        sigma_z2=contr_axi(r_values1, a1, p01)
        print(p0,a)
   
        
        # Tracer les contraintes pour chaque rayon avec une étiquette unique
        label = f'R = 2.5 mm, F = 50 mN, E_polymère = 43 kPA, Sphère-peau_avant-bras'  # Label en mm pour R
        label2 = f'R = 2.5 mm, F = 50 mN, E_polymère = 9 kPA, Sphère-peau_joue'  # Label en mm pour R
        ax1.plot(sigma_z / p0,r_values/a , label=f'σ_z ({label})', color="red")
        ax2.plot(sigma_z2 / p01 ,r_values1/a1 , label=f'σ_z ({label2})', color="orange")
        # ax.plot(-r_values / a, sigma_z / p0, label='_nolegend_', color=color)  # Courbe symétrique
      

    # Configuration des labels et du titre
    ax1.set_title('Contraintes axiales en fonction de la profondeur pour un contact sphère-peau_avant-bras et F = 50mN')
    ax1.set_xlabel('σ_z/p0')
    ax1.set_ylabel('z/a')
    ax1.legend()
    ax1.grid(True)

    # Inverser l'axe des abscisses et le placer à gauche
    ax1.invert_yaxis()
    ax1.invert_xaxis()

    # Inverser les axes
    ax1.xaxis.set_ticks_position('top')
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_label_position('top')
    ax1.yaxis.set_label_position('left')

 


    # Configuration des labels et du titre
    ax2.set_title('Contraintes axiales en fonction de la profondeur pour un contact sphère-peau_joue')
    ax2.set_xlabel('σ_z/p0')
    ax2.set_ylabel('z/a')
    ax2.legend()
    ax2.grid(True)

    # Inverser l'axe des abscisses et le placer à gauche
    ax2.invert_yaxis()
    ax2.invert_xaxis()

    # Inverser les axes
    ax2.xaxis.set_ticks_position('top')
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_label_position('top')
    ax2.yaxis.set_label_position('left')
    
    # Définir les mêmes limites pour l'axe y (z/a) des deux subplots
    y_min = min(ax1.get_ylim()[0], ax2.get_ylim()[0])
    y_max = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(y_min, y_max)
    ax2.set_ylim(y_min, y_max)


    # Afficher le graphique dans le frame Tkinter
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

#####





# Création de la fenêtre principale
root = tk.Tk()
root.title("Analyse des Contraintes")

# Création du widget Notebook pour les onglets
notebook = ttk.Notebook(root)

# Onglet 1 
# tab1 = ttk.Frame(notebook)
# notebook.add(tab1, text='[1]σ_r, R, sp-sp')
# plot_rayon_r_spsp(F_values[0], tab1)

# Onglet 2
# tab2 = ttk.Frame(notebook)
# notebook.add(tab2, text='[2] σ_r, R, sp-pl')
# plot_rayon_r_sppl(F_values[0], tab2)

# Onglet 3
# tab3 = ttk.Frame(notebook)
# notebook.add(tab3, text='[3] σ_θ, R, sp-sp')
# plot_rayon_theta_spsp(F_values[0], tab3)

# Onglet 4 
# tab4 = ttk.Frame(notebook)
# notebook.add(tab4, text='[4] σ_θ, R, sp-pl')
# plot_rayon_theta_sppl(F_values[0], tab4)

# Onglet 5
# tab5 = ttk.Frame(notebook)
# notebook.add(tab5, text='[5] σ_z, R, sp-sp')
# plot_rayon_z_spsp(F_values[0], tab5)

# Onglet 6
tab6 = ttk.Frame(notebook)
notebook.add(tab6, text='[6] σ_z, R, sp-pl')
plot_rayon_z_sppl(F_values[0], tab6)

# Onglet 7
# tab7 = ttk.Frame(notebook)
# notebook.add(tab7, text='[7] σ_r, F, sp-sp')
# plot_force_r_spsp(radii[0], tab7)

# Onglet 8
# tab8 = ttk.Frame(notebook)
# notebook.add(tab8, text='[8] σ_r, F, sp-pl')
# plot_force_r_sppl(radii[0], tab8)

# Onglet 9
# tab9 = ttk.Frame(notebook)
# notebook.add(tab9, text='[9] σ_θ, F, sp-sp')
# plot_force_theta_spsp(radii[0], tab9)

# Onglet 10
# tab10 = ttk.Frame(notebook)
# notebook.add(tab10, text='[10] σ_θ, F, sp-pl')
# plot_force_theta_sppl(radii[0], tab10)

# Onglet 11
# tab11 = ttk.Frame(notebook)
# notebook.add(tab11, text='[11] σ_z, F, sp-sp')
# plot_force_z_spsp(radii[0], tab11)

# Onglet 12
# tab12 = ttk.Frame(notebook)
# notebook.add(tab12, text='[12] σ_z, F, sp-pl')
# plot_force_z_sppl(radii[0], tab12)

# # Onglet 13
# tab13 = ttk.Frame(notebook)
# notebook.add(tab13, text='[12] σ_z, F, sp-pl')
# plot_contrainte_cisaillement(F_values[0], tab13)

# Bouton Quitter
quit_button = tk.Button(root, text='Quitter', command=root.quit)
quit_button.pack(side=tk.BOTTOM, pady=10)


# Affichage des onglets
notebook.pack(fill=tk.BOTH, expand=True)

# Lancement de l'application
root.mainloop()
