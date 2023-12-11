import numpy as np
import numpy.random as rd
from datetime import datetime
import os
import moviepy.video.io.ImageSequenceClip
import argparse
import time
from time import time
import random
# Pour calculer rapidement les distances
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

#==============================================================================
def make_neutron(n, positions, vitesse, sigma):
    """
    Parameters
    ----------
    n : int
        Nombre de nouveaux neutrons par positions.
    positions : numpy.ndarray
        Positions d'où partent les nouveaux neutrons. La matrice doit avoir la
        forme suivante : np.array([[x0, y0, z0], ..., [xi, yi, zi]]).
    vitesse : float
        Vitesse total des neutrons.
    sigma : float
        Ecart-type de la gaussienne qui ajoutera de la variabilité aux
        vitesses des neutrons. Experimetalement < 0,5.

    Returns
    -------
    new_neutrons : numpy.ndarray
        Array 2d contenant les informations des nouveaux neutrons.
        
    """
    n_posi = positions.shape[0]
    # structure des colonnes : t,  x,  y,  z, vx, vy, vz, vtotal
    new_neutrons = np.zeros((n_posi, n, 8))
    # positions de départ
    new_neutrons[:, :, 1:4] +=  positions[:, np.newaxis]
    # vitesse initial
    v_tot = vitesse+np.random.normal(0, sigma, (n_posi, n)) # v total
    rand_colat = np.random.uniform(0, np.pi, (n_posi, n)) # colatitude
    rand_longi = np.random.uniform(0, 2*np.pi, (n_posi, n)) # longitude
    new_neutrons[:, :, 4] = v_tot*np.sin(rand_colat)*np.cos(rand_longi) # vx
    new_neutrons[:, :, 5] = v_tot*np.sin(rand_colat)*np.sin(rand_longi) # vy
    new_neutrons[:, :, 6] = v_tot*np.cos(rand_colat) # vz
    new_neutrons[:, :, 7] = v_tot # v total
    new_neutrons = np.reshape(new_neutrons, (n_posi*n, 8))
    return new_neutrons

def desintegration(data, box_limits, age_fermi=2):
    """
    Parameters
    ----------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.
    box_limits : list
        Limites de la boîte où se situe la pile nucléaire. Contient les
        valeurs caractéristique d'un cylindre : [rayon, hauteur]. La posistion
        du centre du cylindre sera considéré comme étant : (0, 0, h/2). Avec h
        la hauteur du cylindre.
    age_fermi : float, optional
        Age de vie maximale des neutrons. La valeur par deffaut sera 2.

    Returns
    -------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.

    """
    # Pour enlever les neutrons étant hors de la boîte
    data = data[(data[:, 3] >= 0)&(data[:, 3] <= box_limits[1])] # hauteur
    data = data[cdist(np.array([[0, 0]]), data[:, 1:3])[0] <= box_limits[0]] # rayon

    # limite sur la densité local des neutrons
    dmin = 3
    maximum = 30
    distances = cdist(data[:, 1:4], data[:, 1:4])
    for i in range(len(data)):
        nearby = np.argwhere((distances[i] <= dmin)&(
                              distances[i] > 0))[:, 0]

        if nearby.shape[0] > maximum:
            data[nearby[(maximum-1):]] = np.inf

    data = data[data[:, 0] < age_fermi]
    return data

def deplacement(data, delta_t):
    """
    Parameters
    ----------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.
    delta_t : float
        Pas de temps de la simulation.

    Returns
    -------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.

    """
    data[:, 1:4] = data[:, 1:4] + data[:, 4:7] * delta_t
    data[:, 0] = data[:, 0]+delta_t
    return data

def moderation(data, delta_t):
    """
    Parameters
    ----------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.
    delta_t : float
        Pas de temps de la simulation.

    Returns
    -------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.

    """
    mod = 0.9
    # mod ** delta_t => plus dt faible => plus proche de 1
    data[:, 4:] *= (mod ** delta_t) # ||v||*a ~ ||v * a|| avec a très proche de 1
    return data

def fission(data, combustible, min_ray):
    """
    Parameters
    ----------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.
    combustible : numpy.ndarray
        Array 3d listant les positions des particules du combustible.
    min_ray : float
        Distance neutron-noyau en-dessous de laquelle on considère qu'il y a
        une collision (avec ou sans fission).

    Returns
    -------
    data : numpy.ndarray
        Array 2d contenant les informations des neutrons en déplacement.

    """
    n_neuts = data.shape[0]
    n_combs = combustible.shape[0]
    E_th = 1.8 # kj
    distances = cdist(data[:, 1:4], combustible)
    # distance minimal pour chaque neutron par rapport aux particules de
    # combustible. Chaque valeur du vecteur est la distance minimal entre
    # le i-eme neutron et toutes les particules de combustible.
    min_dist = np.argmin(distances, axis=1)
    collider = distances[np.arange(n_neuts), min_dist] <= min_ray
    if len(collider[collider]) > 0:
        # vitesse total des neutrons en collision
        Ec = data[:, 7] # = (vx^2 +vy^2 +vz^2)^0.5
        P = (Ec - E_th) / (Ec * (1 + Ec / (2 * E_th)))

        # True => fission \ False => rebond 
        #  Tirage aléatoire + séléction des neutron 'impactants'
        fission_bool = (np.abs(np.random.normal(0, 1, len(collider))) < P)&(collider)
        rebond_bool = (fission_bool == False)&(collider)

        if len(fission_bool[rebond_bool]) > 0:
            # rebond élastique : m_neutron << m_combustible => combustible ne bouge pas.
            # Pour le moment faire simple :
            data[rebond_bool, 4:7] *= -1

        if len(fission_bool[fission_bool]) > 0: # fission
            new_data = make_neutron(3, data[fission_bool][:, 1:4], 7, 0.2)
            data = np.concatenate((data, new_data))

    return data

def make_combustible(n_couche, delta_ray, pile_ray, pile_height,
                     n_particules, profil='gaussian'):
    """
    Parameters
    ----------
    n_couche : int
        Nombre de couche concentrique de piles.
    delta_ray : float
        Distance radial séparant chaque couche de piles.
    pile_ray : float
        Rayon des piles de combustible.
    pile_height : float
        Hauteur des piles de combustible.
    n_particules : int
        Nombre de particules par pile.
    profil : str, optional
        Type de profil de concentration radiale. La valeur par deffaut
        est 'gaussian'. Les différentes valeures possibles sont : [
        'uniform', 'sqrt', 'gaussian'].

    Returns
    -------
    combustible : numpy.ndarray
        Array 3d listant les positions des particules du combustible.

    """
    # On initialise avec la pile centrale en créant son "centre"
    centre = np.zeros((1, 3))
    # On ajoute les n couches de barres suplémentaires
    for i in range(int(n_couche-1)):
        angles = np.deg2rad(np.arange(0, 360, 60/(i+1)))
        x = np.cos(angles)*(delta_ray*(i+1))
        y = np.sin(angles)*(delta_ray*(i+1))
        z = np.zeros((i+1)*6)
        centre = np.concatenate((centre, np.array([x, y, z]).T))

    if profil == 'uniform':
        rayon_rand = (np.random.uniform(0, 1,
                            (centre.shape[0], n_particules)) **.5)*pile_ray

        theta_rand = np.random.uniform(0, 2*np.pi,
                            (centre.shape[0], n_particules))

        x_rand = np.cos(theta_rand)*rayon_rand
        y_rand = np.sin(theta_rand)*rayon_rand

    elif profil == 'sqrt':
        rayon_rand = (np.random.uniform(0, pile_ray,
                            (centre.shape[0], n_particules)))

        theta_rand = np.random.uniform(0, 2*np.pi,
                            (centre.shape[0], n_particules))

        x_rand = np.cos(theta_rand)*rayon_rand
        y_rand = np.sin(theta_rand)*rayon_rand

    elif profil == 'gaussian':
        x_rand = np.random.normal(0, 1,
                            (centre.shape[0], n_particules))

        y_rand = np.random.normal(0, 1,
                            (centre.shape[0], n_particules))

        theta = np.arctan2(y_rand, x_rand)
        r_rand = (x_rand**2 +y_rand**2)**.5
        r_rand = r_rand/(r_rand[np.arange(centre.shape[0]),
                                np.argmax(r_rand, axis=1),
                                np.newaxis])*pile_ray

        x_rand = r_rand*np.cos(theta)
        y_rand = r_rand*np.sin(theta)
        

    z_rand = np.random.uniform(0, pile_height,(centre.shape[0], n_particules))

    p_rand = np.array([x_rand, y_rand, z_rand]).transpose((1, 2, 0))

    combustible = p_rand+centre[:, np.newaxis]
    combustible = np.concatenate(combustible)
    return combustible

def movie(i, dossier):
    image_folder = "./"+dossier+"/"
    fps = 4
    image_files = []
    for n in range(1, i+1):
        image_files.append(os.path.join(image_folder, str(n)+'p.png'))

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(image_folder+'ball-p.mp4')
    for n in range(1, i+1):
        os.remove('./'+dossier+"/"+str(n)+"p.png")

def cycle(n_steps, dt, min_ray, combustible, neutrons, box_limits, dossier, scale, width, height):
    start = datetime.now()
    kmoy = 0
    keff = []
    neutron = []
    neutron.append(neutrons.shape[0])
    projected_comb = project_points(combustible[:, 0], combustible[:, 1], scale, width, height)
    for i in range(1, int(n_steps+1)):
        neutrons = deplacement(neutrons, dt)
        neutrons = fission(neutrons, combustible, min_ray)
        neutrons = moderation(neutrons, dt)
        neutrons = desintegration(neutrons, box_limits)
        length = neutrons.shape[0]
        if length > 0:
            neutron.append(length)
            k = neutron[i]/neutron[i-1]
            kmoy = kmoy + k
            keff.append(kmoy/i)
            delta = datetime.now() - start
            print("Cycle n°"+str(i))
            print("Il y a "+str(neutron[i])+" neutrons")
            print("k = "+str(k))
            print("kmoy = "+str(keff[i-1]))
            print ("Temps d'exécution : ", delta, '\n')

            projected_points = project_points(neutrons[:, 1], neutrons[:, 2], scale, width, height)

            image = Image.new("RGB", (width, height), (0, 0, 0))
            draw = ImageDraw.Draw(image)
            draw.point(tuple(map(tuple, projected_points)), fill=(0, 255, 255))
            draw.point(tuple(map(tuple, projected_comb)), fill=(0, 255, 0))
            image.save('./'+dossier+"/"+str(i)+"p.png")

        else:
            break

    delta = datetime.now() - start
    print ("Temps d'exécution : ", delta)
    movie(i, dossier)
    fig = plt.figure()
    plt.plot(range(len(neutron)), neutron)
    plt.title("Variation de Neutrons", fontsize=14)
    plt.savefig('./'+dossier+'/var-p.png')
    plt.close('all')

    fig = plt.figure()
    plt.plot(range(len(keff)), keff)
    plt.title("Variation de keff", fontsize=14)
    plt.savefig('./'+dossier+'/keff-p.png')
    plt.close('all')

def compet(n_steps, dt, delta_ray, n_couche, pile_ray, pile_height, n_particules, min_ray, scale, width, height, dossier):
    if not os.path.exists(dossier):
        os.mkdir(dossier)

    combustible = make_combustible(n_couche, delta_ray, pile_ray, pile_height,n_particules)

    data = make_neutron(100, np.array([[0, 0, pile_height/2]]), 7, 0.2)
    box_limits = [(n_couche*delta_ray + pile_ray/2), pile_height]
    cycle(n_steps, dt, min_ray, combustible, data, box_limits, dossier, scale, width, height)

def project_points(x, y, scale, width, height):
    proj = np.column_stack((x * scale + width/2, y * scale + height/2))
    return proj

width, height = 800, 800
scale = 10.

parser = argparse.ArgumentParser(description='PHP-DEM')
parser.add_argument('n', help='ex 200')
parser.add_argument('dt', help='ex 0.1')
parser.add_argument('dr', help='ex 10')
parser.add_argument('rep', help='ex 3')
parser.add_argument('min_r', help='ex 1')
parser.add_argument('pile_h', help='ex 12')
parser.add_argument('pile_r', help='ex 6')
parser.add_argument('n_comb', help='ex 2')
parser.add_argument('dossier', help='ex toto')
args = parser.parse_args()

compet(int(args.n), float(args.dt), float(args.dr), float(args.rep), float(args.pile_r), float(args.pile_h), int(args.n_comb)*100, float(args.min_r), scale, width, height, str(args.dossier))
