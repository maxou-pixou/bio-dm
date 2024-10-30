import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Paramètres physiques
F_values = [50, 100]  # Forces appliquées (en Newtons)
radii = [10e-3, 5e-3, 2e-3, 500e-6]  # Différents rayons de la sphère (en mètres)
E_sphere = 210e9  # Module de Young de la sphère (acier, en Pascals)
nu_sphere = 0.3  # Coefficient de Poisson de la sphère
E_polymer_high = 10e6  # Module de Young du polymère (en Pascals)
nu_polymer = 0.45  # Coefficient de Poisson du polymère
E_values =  [10e6, 210e9]


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
    sigma_r1 = (1 - 2 * nu / 3) * ((a / r)**2) * \
              (1 - (1 - (r / a)**2)**(3/2)) - (1 - (r / a)**2)**(1/2)
    sigma_r = sigma_r1 / p0

    sigma_theta1 = -(1 - 2 * nu / 3) * ((a / r)**2) * \
                  (1 - (1 - (r / a)**2)**(3/2)) - 2 * nu * (1 - (r / a)**2)**(1/2)
    sigma_theta = sigma_theta1 / p0
    
    sigma_z1 = - (1-(r / a)**(2))**1/2
    sigma_z = sigma_z1 / p0

    return sigma_r, sigma_theta, sigma_z

# Fonction pour calculer les contraintes radiales, circonférentielles et axiales EXTERIEURES
def contraintes_polairesext(r, a, p0, nu):
    sigma_r2 = (1 - 2 * nu) * (a **2) * \
              3*r**2
    sigma_r3 = sigma_r2 / p0

    sigma_theta3 = -sigma_r3 / p0

    return sigma_r3, sigma_theta3

# Fonction pour tracer les contraintes en fonction du rayon
def plot_constraints_vs_radius(F, frame):
    global r2
    E_star = reduced_modulus(E_sphere, nu_sphere, E_polymer_high, nu_polymer)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    num_radii = len(radii)
    colors = plt.cm.viridis(np.linspace(0, 1, num_radii))  # Utilisation de la colormap viridis

    for i, R in enumerate(radii):
        a = contact_radius(F, R, E_star)
        
        r = np.linspace(0.01*a,a, 100)
        p0 = max_pressure(F, a)
        
        list_sr = []
        list_sz = []

        sigma_r, sigma_theta, sigma_z = contraintes_polaires(r, a, p0, nu_polymer)

        list_sr.append(sigma_r)
        list_sz.append(sigma_z)
       
        # Calcul du score R² pour la régression linéaire de sigma_r
        X = r.reshape(-1, 1)  # Reshape pour sklearn
        y1 = sigma_r / 1e6
        y2 = sigma_theta / 1e6
        y3 = sigma_z / 1e6

        model1 = LinearRegression().fit(X, y1)
        model2 = LinearRegression().fit(X, y2)
        model3 = LinearRegression().fit(X, y3)

        y_pred1 = model1.predict(X)
        y_pred2 = model2.predict(X)
        y_pred3 = model3.predict(X)

        r2_r = round(r2_score(y1, y_pred1),2)
        r2_theta = round(r2_score(y2, y_pred2),2)
        r2_z = round(r2_score(y3, y_pred3),2)

        #label de chaque courbe
        label = f'R = {R} m, F = {F} N'
        # Tracé des contraintes et des contraintes symétrisées avec la même couleur
        color = colors[i]
        line_r, = ax.plot((r)/a, sigma_r / 1e6, label=f'σ_r ({label}), R² = {r2_r}', color=color)
        ax.plot(-(r)/a, sigma_r / 1e6, '--', color=color, label='_nolegend_')  # σ_r en MPa
        
        line_theta, = ax.plot((r)/a, sigma_theta / 1e6, label=f'σ_θ ({label}), R² = {r2_theta}', color=color)
        ax.plot(-(r)/a, sigma_theta / 1e6, '--', color=color, label='_nolegend_')  # σ_θ en MPa
        
        line_z, = ax.plot((r)/a, sigma_z / 1e6, label=f'σ_z ({label}), R² = {r2_z}', color=color)
        ax.plot(-(r)/a, sigma_z / 1e6, '--', color=color, label='_nolegend_')  # σ_z en MPa

        #Ajouter des annotations avec des chiffres à côté des courbes
        ax.annotate(f'σ_r, R={R} m', xy=((r[-1])/a, sigma_r[-1] / 1e6), color=color)
        ax.annotate(f'σ_θ, R={R} m', xy=((r[-1])/a, sigma_theta[-1] / 1e6), color=color)
        ax.annotate(f'σ_z, R={R} m', xy=((r[-1])/a, sigma_z[-1] / 1e6), color=color)

        

    ax.set_title('Contraintes intérieures pour différents rayons pour F = 50 N')
    ax.set_xlabel('Rayon de la zone de contact (mm)')
    ax.set_ylabel('Contraintes (MPa)')
    ax.legend()
    ax.grid(True)
  
    ax.invert_yaxis()

    # Afficher le graphique dans l'onglet
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    return list_sz, list_sr
# Fonction pour tracer les contraintes en fonction de la force
def plot_constraints_vs_force(r, E_star, frame):
    fig, ax = plt.subplots(figsize=(8, 6))
    num_forces = len(F_values)
    colors = plt.cm.viridis(np.linspace(0, 1, num_forces))  # Utilisation de la colormap viridis
    
    

    for i, F in enumerate(F_values):
        a = contact_radius(F, r, E_star)
        r_vals = np.linspace(0.01 * a, a, 100)
        p0 = max_pressure(F, a)
        sigma_r, sigma_theta, sigma_z = contraintes_polaires(r_vals, a, p0, nu_polymer)

        list_sr2 = []
        list_sz2 = []

        list_sr2.append(sigma_r)
        list_sz2.append(sigma_z)
        
      
        
        # Calcul du score R² pour la régression linéaire de sigma_r
        X = r_vals.reshape(-1, 1)  # Reshape pour sklearn
        y1 = sigma_r / 1e6
        y2 = sigma_theta / 1e6
        y3 = sigma_z / 1e6

        model1 = LinearRegression().fit(X, y1)
        model2 = LinearRegression().fit(X, y2)
        model3 = LinearRegression().fit(X, y3)

        y_pred1 = model1.predict(X)
        y_pred2 = model2.predict(X)
        y_pred3 = model3.predict(X)

        r2_r = round(r2_score(y1, y_pred1),2)
        r2_theta = round(r2_score(y2, y_pred2),2)
        r2_z = round(r2_score(y3, y_pred3),2)
       

        #label de chaque courbe
        label = f'F = {F} N, R = {10} mm'
        
        # Tracé des contraintes et des contraintes symétrisées avec la même couleur
        color = colors[i]
        line_r, = ax.plot((r_vals)/a, sigma_r / 1e6, label=f'σ_r ({label}), R²={r2_r}', color=color)
        ax.plot(-(r_vals)/a, sigma_r / 1e6, '--', color=color, label='_nolegend_')  # σ_r en MPa
        
        line_theta, = ax.plot((r_vals)/a, sigma_theta / 1e6, label=f'σ_θ ({label}), R²={r2_theta}', color=color)
        ax.plot(-(r_vals)/a, sigma_theta / 1e6, '--', color=color, label='_nolegend_')  # σ_θ en MPa
        
        line_z, = ax.plot((r_vals)/a, sigma_z / 1e6, label=f'σ_z ({label}), R²={r2_z}', color=color)
        ax.plot(-(r_vals)/a, sigma_z / 1e6, '--', color=color, label='_nolegend_')  # σ_z en MPa

        # Ajouter des annotations avec la contrainte polaire et la valeur de la force à côté des courbes
        ax.annotate(f'σ_r, F={F} N', xy=((r_vals[-1])/a, sigma_r[-1] / 1e6), color=color)
        ax.annotate(f'σ_θ, F={F} N', xy=((r_vals[-1])/a, sigma_theta[-1] / 1e6), color=color)
        ax.annotate(f'σ_z, F={F} N', xy=((r_vals[-1])/a, sigma_z[-1] / 1e6), color=color)

    ax.set_title('Contraintes intérieures en fonction du rayon de contact pour une bille de rayon R = 10 mm')
    ax.set_xlabel('Rayon de la zone de contact (mm)')
    ax.set_ylabel('Contraintes (MPa)')
    ax.legend()
    ax.grid(True)
    ax.invert_yaxis()

    # Afficher le graphique dans l'onglet
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    return list_sz2, list_sr2

# Fonction pour tracer les contraintes en fonction du rayon
def plot_constraints_vs_radiusext(F, frame):
    E_star = reduced_modulus(E_sphere, nu_sphere, E_polymer_high, nu_polymer)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    num_radii = len(radii)
    colors = plt.cm.viridis(np.linspace(0, 1, num_radii))  # Utilisation de la colormap viridis

    for i, R in enumerate(radii):
        a = contact_radius(F, R, E_star)
        r = np.linspace(0.01 * a, a, 100)
        p0 = max_pressure(F, a)
        sigma_r3, sigma_theta3 = contraintes_polairesext(r, a, p0, nu_polymer)

         # Calcul du score R² pour la régression linéaire de sigma_r
        X = r.reshape(-1, 1)  # Reshape pour sklearn
        y1 = sigma_r3 / 1e6
        y2 = sigma_theta3 / 1e6
       
        model1 = LinearRegression().fit(X, y1)
        model2 = LinearRegression().fit(X, y2)
     
        y_pred1 = model1.predict(X)
        y_pred2 = model2.predict(X)
       
        r2_r3 = round(r2_score(y1, y_pred1),2)
        r2_theta3 = round(r2_score(y2, y_pred2),2)
      
        #label de chaque courbe
        label = f'R = {R*1e3} mm, F = {F} N'
        
        # Tracé des contraintes et des contraintes symétrisées avec la même couleur
        color = colors[i]
        line_r, = ax.plot((r)/a, sigma_r3 / 1e6, label=f'σ_r ({label}), R² = {r2_r3}', color=color)
        ax.plot(-(r)/a, sigma_r3 / 1e6, '--', color=color, label='_nolegend_')  # σ_r en MPa
        
        line_theta, = ax.plot((r)/a, sigma_theta3 / 1e6, label=f'σ_θ ({label}), R² = {r2_theta3}', color=color)
        ax.plot(-(r)/a, sigma_theta3 / 1e6, '--', color=color, label='_nolegend_')  # σ_θ en MPa

        #Ajouter des annotations avec des chiffres à côté des courbes
        ax.annotate(f'σ_r, R={R} m', xy=((r[-1])/a, sigma_r3[-1] / 1e6), color=color)
        ax.annotate(f'σ_θ, R={R} m', xy=((r[-1])/a, sigma_theta3[-1] / 1e6), color=color)
       

    ax.set_title('Contraintes intérieures en fonction du rayon de contact pour F = 50 N')
    ax.set_xlabel('Rayon de la zone de contact (mm)')
    ax.set_ylabel('Contraintes (MPa)')
    ax.legend()
    ax.grid(True)



    # Afficher le graphique dans l'onglet
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


# Fonction pour tracer les contraintes en fonction de la force
def plot_constraints_vs_forceext(r, E, frame):
    
    fig, ax = plt.subplots(figsize=(8, 6))
    num_forces = len(F_values)
    colors = plt.cm.viridis(np.linspace(0, 1, num_forces))  # Utilisation de la colormap viridis

    for i, F in enumerate(F_values):
        a = contact_radius(F, r, E)
        r_vals = np.linspace(0.01 * a, a, 100)
        p0 = max_pressure(F, a)
        sigma_r3, sigma_theta3= contraintes_polairesext(r_vals, a, p0, nu_polymer)

        # Calcul du score R² pour la régression linéaire de sigma_r
        X = r_vals.reshape(-1, 1)  # Reshape pour sklearn
        y1 = sigma_r3 / 1e6
        y2 = sigma_theta3 / 1e6
       
        model1 = LinearRegression().fit(X, y1)
        model2 = LinearRegression().fit(X, y2)

        y_pred1 = model1.predict(X)
        y_pred2 = model2.predict(X)

        r2_r3 = round(r2_score(y1, y_pred1),2)
        r2_theta3 = round(r2_score(y2, y_pred2),2)
      
       #label des courbes
        label = f'F = {F} N, r = {10} mm'
        
        # Tracé des contraintes et des contraintes symétrisées avec la même couleur
        color = colors[i]
        line_r, = ax.plot((r_vals)/a, sigma_r3 / 1e6, label=f'σ_r ({label}), R² = {r2_r3}', color=color)
        ax.plot(-(r_vals)/a, sigma_r3 / 1e6, '--', color=color, label='_nolegend_')  # σ_r en MPa
        
        line_theta, = ax.plot((r_vals)/a, sigma_theta3 / 1e6, label=f'σ_θ ({label}), R² = {r2_theta3}', color=color)
        ax.plot(-(r_vals)/a, sigma_theta3 / 1e6, '--', color=color, label='_nolegend_')  # σ_θ en MPa
        
        # Ajouter des annotations avec la contrainte polaire et la valeur de la force à côté des courbes
        ax.annotate(f'σ_r, F={F} N', xy=((r_vals[-1])/a, sigma_r3[-1] / 1e6), color=color)
        ax.annotate(f'σ_θ, F={F} N', xy=((r_vals[-1])/a, sigma_theta3[-1] / 1e6), color=color)
        
    ax.set_title('Contraintes intérieures en fonction de la force pour R = 1 m, E = 10 MPa')
    ax.set_xlabel('Rayon de la zone de contact (mm)')
    ax.set_ylabel('Contraintes (MPa)')
    ax.legend()
    ax.grid(True)

    # Afficher le graphique dans l'onglet
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

def plot_cisaillement_radius(f,frame):
    fig, ax = plt.subplots(figsize=(8, 6))
    E_star = reduced_modulus(E_sphere, nu_sphere, E_polymer_high, nu_polymer)
   
    num_radii=len(radii)
    colors = plt.cm.viridis(np.linspace(0, 1, num_radii))
    for i, R in enumerate(radii):
        a = contact_radius(f, R, E_star)
        p = max_pressure(f,a)
        r = np.linspace(0.01*a,a,100)
        z = 0.48*a
        sigma_r4,sigma_theta4,sigma_z4=contraintes_polaires(r, a, p, nu_polymer)
        
        z_lin=np.linspace(0.01*z,z,100)

        #tau_calcul=0.31*p
        tau_calcul=(sigma_z4-sigma_r4)/2
        
        print(z_lin)
        
        color=colors[i]   
        ax.plot(r/a, z_lin/a, label=f'R = {R} mm, τ = {tau_calcul} MPa', color=color)
    
    ax.set_title('Contraintes de cisaillement en fonction du rayon')
    ax.set_xlabel('Rayon de la zone de contact (r/a)')
    ax.set_ylabel('Contraintes de cisaillement (z/a)')
    ax.legend()
    ax.grid(True)

    # Afficher le graphique dans l'onglet
    canvas = FigureCanvasTkAgg(fig, master=frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    

# Création de la fenêtre principale
root = tk.Tk()
root.title("Analyse des Contraintes")

# Création du widget Notebook pour les onglets
notebook = ttk.Notebook(root)

# Onglet 1 : Contraintes int en fonction du rayon
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text='Rayon (int)')
plot_constraints_vs_radius(F_values[0], tab1)

# Onglet 2 : Contraintes int en fonction de la force
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text='Force (int)')
plot_constraints_vs_force(1.0, E_polymer_high, tab2)  



# Onglet 3 : Contraintes en ext fonction du module de Young
tab3 = ttk.Frame(notebook)
notebook.add(tab3, text='Rayon (ext)')
plot_constraints_vs_radiusext(F_values[0],  tab3)  

# Onglet 4 : Contraintes en ext fonction du module de Young
tab4= ttk.Frame(notebook)
notebook.add(tab4, text='Force (ext)')
plot_constraints_vs_forceext(1.0,E_polymer_high,  tab4) 

# Onglet 5 : Contraintes de cisaillement
tab5 = ttk.Frame(notebook)
notebook.add(tab5, text='Contrainte de Cisaillement pour différents rayon à force F = 50N')
plot_cisaillement_radius(F_values[0], tab5)

# Bouton Quitter
quit_button = tk.Button(root, text='Quitter', command=root.quit)
quit_button.pack(side=tk.BOTTOM, pady=10)


# Affichage des onglets
notebook.pack(fill=tk.BOTH, expand=True)

# Lancement de l'application
root.mainloop()
