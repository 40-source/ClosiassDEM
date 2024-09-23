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

def histoire(position,direction,e):
    centre = np.array([0, 0, 0]) #On fait l'hypothèse d'un coeur homogène
    rayon_echantillonne = 153 # spécifique rayon d'un coeur

    if(distance(centre,position,rayon_echantillonne)):
        N_comb = N(18.7,235)*5/100 # 235 à 5%

        Ss = section_ef(N_comb,section_s(e))
        Sc = section_ef(N_comb,section_c(e))
        Sf = section_ef(N_comb,section_f(e))

        S = Ss + Sc + Sf
        
        chi = np.random.uniform(0,1)
        longueur = s(S,chi) # cm
        position = position.copy() + np.array(direction) * longueur

        if(chi<Ss/S):
            return diffusion(position,direction,e,235),[] #diffusion du noyau U235
        elif(chi>(Ss+Sc)/S):
            return [],[np.round(position,0)] #fission du noyau
        else:
            return [],[]#capture stérile
    else:
        S = 0.931*10**-27*N(1,1) #Eau légère
        chi = np.random.uniform(0,1)
        longueur = s(S,chi) # cm
        position = position.copy() + np.array(direction) * longueur

        return diffusion(position,direction,e,1),[] #diffusion H1

def distance(centre, position, rayon_echantillonne):
    squared_distance = np.sum((np.array(position) - np.array(centre))**2)
    return squared_distance < rayon_echantillonne**2

def present(position,ensemble):
    if len(position) > 0:
        for i in ensemble:
            if len(i)>0:
                if i[0][0] == position[0][0] and i[0][1] == position[0][1] and i[0][2] == position[0][2]:
                    return True
        return False
    else:
        return False

def bach(source, n):
    rayon_echantillonne=180 #générale
    centre = np.array([0, 0, 0])
    red = []
    purple = []
    for i in range(n):
        future = []
        for neutron in source:
            if distance(centre, neutron[0], rayon_echantillonne):
                h = histoire(neutron[0], neutron[1], neutron[2])
                future.append(h[0])  
                if present(h[1],red):
                    purple.append(h[1])
                else:
                    red.append(h[1])
        future = [item for sublist in future for item in sublist]
        red_f = [item for sublist in red for item in sublist]
        purple_f = [item for sublist in purple for item in sublist]

        source = future.copy()
        for _ in range(200):
            source.append([[0,0,0],direction(),2])

        if len(future) > 0:
            projected_points = project_points(np.array([row[0] for row in future]))
        else:
            projected_points = []
        
        if len(red) > 0:
            project_points_red = project_points(np.array(red_f))
        else:
            project_points_red = []
        
        if len(purple) > 0:
            project_points_purple = project_points(np.array(purple_f))
        else:
            project_points_purple = []

        image = create_image(projected_points,project_points_red,project_points_purple)
        image.save(f"{i}.png")

        print(len(future),len(red),len(purple))

    movie(i)

def create_image(projected_points,red,purple):
    image = Image.new("RGB", (width, height), (0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.point(tuple(map(tuple, projected_points)), fill=(0, 255, 255))
    draw.point(tuple(map(tuple, red)), fill=(255, 0, 0))
    draw.point(tuple(map(tuple, purple)), fill=(255, 255, 255))
    return image

def project_points(points):
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    return np.column_stack((x * scale + width/2, y * scale + height/2))

def movie(i):
    image_folder = "./"
    fps = 4
    image_files = [os.path.join(image_folder, f"{n}.png") for n in range(i + 1)]

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(os.path.join(image_folder, 'MC.mp4'))

    for file in image_files:
        os.remove(file)

width, height = 2400, 2400
scale = 1

def start(m,n):
    source=[]
    for _ in range(m):
        source.append([[0,0,0],direction(),2])
    bach(source,n)

start(500,120)