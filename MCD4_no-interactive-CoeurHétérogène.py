import numpy as np
import os
from scipy.interpolate import interp1d
from PIL import Image, ImageDraw
import moviepy.video.io.ImageSequenceClip

energies = [0.001,0.025, 1.0, 2.0, 5.0, 10.0]  # MeV

sigma_f_values = np.log([10,3, 1, 1.05, 1.07, 4])  # barns
interp_f = interp1d(energies, sigma_f_values, kind='linear', fill_value='extrapolate')

sigma_c_values = np.log([5,0.5, 0.1, 0.05, 0.001, 0.0001])  # barns
interp_c = interp1d(energies, sigma_c_values, kind='linear', fill_value='extrapolate')

sigma_s_values = np.log([1, 1.01, 1.02, 1.04, 4, 4.23])  # barns
interp_s = interp1d(energies, sigma_s_values, kind='linear', fill_value='extrapolate')

def section_s(e):
    return np.exp(interp_s(e))*10**-24

def section_c(e):
    return np.exp(interp_c(e))*10**-24

def section_f(e):
    return np.exp(interp_f(e))*10**-24

def N(p,A):
    return 6.023*10**23*p/A

def section_ef(N_comb,s):
    return N_comb*s

def s(section,chi):
    return -1/section*np.log(chi)

def direction():
    direction_vector = np.random.uniform(-1, 1, 3)
    norme = np.linalg.norm(direction_vector)
    return direction_vector / norme if norme != 0 else direction_vector

def diffusion(position, direction_p, e, A):
    new_direction = np.array(direction())
    cos_theta = np.dot(new_direction, direction_p)
    e = e * (A**2 + 2*A*cos_theta + 1) / (A + 1)**2
    return [[position, new_direction, e]]

def fission(position):
    resultat = []
    e = 2 # 2 Mev
    n = int(np.random.uniform(2,4)) # espérence proche de 2.4
    for _ in range(n):
        new_direction = direction()
        new_position = position.copy()
        resultat.append([position,new_direction,e])
    return resultat

def carto(cart, position):
    rayon_echantillonne = 20 # spécifique
    for centre in cart:
        if distance_carrer(centre,position,rayon_echantillonne):
            return True
    return False


def histoire(position, direction, e, cart):
    if(carto(cart, position)):
        N_comb = N(18.7,235)*5/100 # 235 à 5%

        Ss = section_ef(N_comb,section_s(e))
        Sc = section_ef(N_comb,section_c(e))
        Sf = section_ef(N_comb,section_f(e))

        S = Ss + Sc + Sf
        
        chi = np.random.uniform(0,1)
        longueur = s(S,chi) # cm
        position = position.copy() + np.array(direction) * longueur

        if(chi<Ss/S):
            return diffusion(position,direction,e,235) #diffusion du noyau U235
        elif(chi>(Ss+Sc)/S):
            return fission(position) #fission du noyau
        else:
            return []#capture stérile
    else:
        S = 0.931*10**-27*N(1,1) #Eau légère
        chi = np.random.uniform(0,1)
        longueur = s(S,chi) # cm
        position = position.copy() + np.array(direction) * longueur

        return diffusion(position,direction,e,1) #diffusion H1

def distance(centre, position, rayon_echantillonne):
    squared_distance = np.sum((np.array(position) - np.array(centre))**2)
    return squared_distance < rayon_echantillonne**2

def distance_carrer(centre, position, rayon_echantillonne):
    i = 0
    hauteur = 100
    return abs(position[0]-centre[0]) < rayon_echantillonne and abs(position[1]-centre[1]) < rayon_echantillonne and abs(position[2]-centre[2]) < hauteur

def bach(source, n, cart):
    rayon_echantillonne=180 #générale
    centre = np.array([0, 0, 0])
    for i in range(n):
        future = [histoire(neutron[0], neutron[1], neutron[2], cart) for neutron in source if distance(centre[:2], neutron[0][:2], rayon_echantillonne)]
        future = [item for sublist in future for item in sublist]

        source = future.copy()

        if len(future) > 0:
            projected_points = project_points(np.array([row[0] for row in future]))
        else:
            projected_points = []

        image = create_image(projected_points)
        image.save(f"{i}.png")

        print(len(future))

    movie(i)

def create_image(projected_points):
    image = Image.new("RGB", (width, height), (0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.point(tuple(map(tuple, projected_points)), fill=(0, 255, 255))
    return image

def project_points(points):
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    return np.column_stack((x * scale + width/2, y * scale + height/2))

def movie(i):
    image_folder = "./"
    fps = 20
    image_files = [os.path.join(image_folder, f"{n}.png") for n in range(i + 1)]

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(os.path.join(image_folder, 'MC.mp4'))

    for file in image_files:
        os.remove(file)

def coeur(l):
    rayon_echantillonne = 20
    indices = np.arange(-l, l, rayon_echantillonne)
    i, j = np.meshgrid(indices, indices)
    mask = np.sqrt(i**2 + j**2) < l
    cart = np.column_stack((i[mask], j[mask], np.zeros(mask.sum())))
    return cart

width, height = 1200, 1200
scale = 2

def start(m,n,cart):
    source=[]
    for _ in range(m):
        source.append([[0,0,0],direction(),2])
    bach(source,n,cart)

cart = coeur(153)
print(cart)
start(20000,620,cart)