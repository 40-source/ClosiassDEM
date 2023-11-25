import numpy as np
import numpy.random as rd
from datetime import datetime
import os
import moviepy.video.io.ImageSequenceClip
import argparse
import time
import random
from numba import njit
import numba as nb
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

#path = "/home/TuanWang/Desktop/Sim/"
#os.chdir(path)

def neutronvect(N, ata0, ata1, ata2, position):
    for i in range(N):
        t = [0.0]
        v = 10
        vitesse=np.array([rd.uniform(-1,1),rd.uniform(-1,1),rd.uniform(-1,1)])
        vitesse = vitesse * v / np.linalg.norm(vitesse)
        ata0 = np.append(ata0, [position], axis=0)
        ata1 = np.append(ata1, [vitesse], axis=0)
        ata2 = np.append(ata2, [t], axis=0)
    return (ata0, ata1, ata2)

@njit
def desintegration_aux(ata0, ata2):
    dmin = 3
    bulle_max = 30
    for i in range(len(ata0)):
        dif = ata0 - ata0[i]
        distances = np.sqrt(np.sum(dif*dif, axis=1))
        nearby_indices = np.where((distances > 0)&(distances <= dmin))[0]
        ata2[nearby_indices[bulle_max:], 0] = np.inf

def desintegration(ata0, ata1, ata2, dinf, d):
    age_fermi = 2
    desintegration_aux(ata0, ata2)
    mask = (np.sqrt(ata0[:,1]**2 + ata0[:,1]**2 + (ata0[:,1]-d/2)**2 ) < dinf) & (np.ravel(ata2 <= age_fermi))
    ata0 = ata0[mask]
    ata1 = ata1[mask]
    ata2 = ata2[mask]
    return (ata0, ata1, ata2)

@njit
def deplacement(ata0, ata1, ata2, t):
    ata0 += ata1 * t
    ata2 += t
    return (ata0,ata2)

@njit
def moderation(ata1, t):
    mod = 0.9
    ata1 *= mod ** t
    return ata1

@njit
def fission_aux(ata1,ata2,i):
    ata1[i]*=rd.normal(-1,1)
    ata1[i]*=rd.normal(-1,1)
    ata1[i]*=rd.normal(-1,1)
    ata2[i][0]=0

def fission(ata0, ata1, ata2, comb, r):
    i=0
    V0=1.8
    m=0.1
    for i in range(len(ata0)):
        distances = np.linalg.norm(comb - ata0[i],axis=1)
        nearby_indices = np.where(distances <= r)[0]
        if nearby_indices.size == 0:
            i += 1
            continue
        noyau = np.argmin(distances)
        Ec=np.sum(ata1[i]**2)*m
        if Ec < V0:
            R=0
        else:
            R=((np.sqrt(Ec)-np.sqrt(Ec)*np.sqrt(Ec-V0))/(np.sqrt(Ec)+np.sqrt(Ec)*np.sqrt(Ec-V0)))**2
        if(rd.normal(0,1)>R):
            ata0, ata1, ata2 = neutronvect(3, ata0, ata1, ata2, comb[noyau])
            ata2[i][0]=np.inf
        else:
            fission_aux(ata1,ata2,i)
    return ata0, ata1, ata2

def combustible(r,d,n,e,rep):
    comb=np.array([[0.0,0.0,d/2]])
    k=1
    while(k<rep):
        i=0
        angle=2*np.pi/(n*k)
        while i<n*k+1:
            x=np.cos(angle*i)*r*k
            y=np.sin(angle*i)*r*k

            if(i==n and k==1):
                x=0
                y=0
            j=0
            while j<d:
                f=1
                while f<e:
                    xf=rd.normal(x,1)
                    yf=rd.normal(y,1)
                    comb = np.append(comb, [[xf,yf,j]], axis=0)
                    f=f+1
                j=j+0.1
            i=i+1
        k=k+1
    return comb

def movie(i,dossier):
    image_folder="./"+dossier+"/"
    fps=4
    image_files=[]
    n=1
    while n<=i:
        image_files.append(os.path.join(image_folder,str(n)+'p.png'))
        n=n+1
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(image_folder+'ball-p.mp4')
    n=1
    while n<=i:
        os.remove('./'+dossier+"/"+str(n)+"p.png")
        n=n+1

def cycle(n,t,r,comb,ata0,ata1,ata2,dinf,d,dossier):
    start = datetime.now()
    i=0
    tmp=len(ata0)
    kmoy=0
    keff=[]
    neutron=[]
    neutron.append(len(ata0))
    projected_comb = project_points(comb)
    while i<n:
        ata0,ata2=deplacement(ata0,ata1,ata2,t)
        ata0,ata1,ata2=fission(ata0,ata1,ata2,comb,r)
        ata1=moderation(ata1,t)
        ata0,ata1,ata2=desintegration(ata0, ata1, ata2, dinf, d)
        i=i+1
        print("Cycle n°"+str(i))
        k=len(ata0)/tmp
        kmoy=kmoy+k
        print("Il y a "+str(len(ata0))+" neutrons")
        neutron.append(len(ata0))
        print("k = "+str(k))
        print("kmoy = "+str(kmoy/i))
        tmp=len(ata0)
        delta=datetime.now() - start
        print ("Temps d'exécution : ", delta)
        keff.append(kmoy/i)
        image = Image.new("RGB", (width, height), (0, 0, 0))
        draw = ImageDraw.Draw(image)
        projected_points = project_points(ata0)
        draw.point(tuple(map(tuple, projected_points)), fill=(0, 255, 255))
        draw.point(tuple(map(tuple, projected_comb)), fill=(0, 255, 0))
        image.save('./'+dossier+"/"+str(i)+"p.png")
    delta=datetime.now() - start
    print ("Temps d'exécution : ", delta)
    movie(i,dossier)
    fig = plt.figure()
    plt.plot([i for i in range(len(neutron))],neutron)
    plt.title("Variation de Neutrons",fontsize = 14)
    plt.savefig('./'+dossier+'/var-p.png')
    plt.close('all')
    fig = plt.figure()
    plt.plot([i for i in range(len(keff))],keff)
    plt.title("Variation de keff",fontsize = 14)
    plt.savefig('./'+dossier+'/keff-p.png')
    plt.close('all')

def compet(n,t,r,rep,uma,d,b,e):
    n=float(n)
    t=float(t)
    r=float(r)
    rep=float(rep)
    uma=float(uma)
    d=float(d)
    b=float(b)
    e=float(e)
    dossier=str(n)+str(t)+str(r)+str(rep)+str(uma)+str(d)+str(b)+str(e)
    if not os.path.exists(dossier):
        os.mkdir(dossier)
    N = 100
    position = [0.0, 0.0, d/2]  # Valeur initiale de position
    ata0 = np.array([position])  # Initialisation de ata0 avec la première position
    ata1 = np.array([[0.0,0.0,0.0]])  # Crée un tableau vide pour ata1
    ata2 = np.array([[0.0]])  # Crée un tableau vide pour ata2

    ata0,ata1,ata2=neutronvect(N, ata0, ata1, ata2, position)
    cycle(n,t,uma,combustible(r,d,b,e,rep),ata0,ata1,ata2,r*rep,d,dossier)

@njit
def project_points(points):
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    return np.column_stack((x * scale + width/2, y * scale + height/2))

width, height = 800, 800
scale = 10
#compet(120,0.1,10,3,1,12,6,2)

parser = argparse.ArgumentParser(description='PHP-DEM')
parser.add_argument('n', help='ex 100')
parser.add_argument('t', help='ex 0.1')
parser.add_argument('r', help='ex 10')
parser.add_argument('rep', help='ex 3')
parser.add_argument('uma', help='ex 1')
parser.add_argument('d', help='ex 12')
parser.add_argument('b', help='ex 6')
parser.add_argument('e', help='ex 2')
args = parser.parse_args()
compet(args.n,args.t,args.r,args.rep,args.uma,args.d,args.b,args.e)