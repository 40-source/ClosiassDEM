import numpy as np
import numpy.random as rd
from datetime import datetime
import os
import moviepy.video.io.ImageSequenceClip
import argparse
import time
import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import random
from numba import njit
import numba as nb
import matplotlib.pyplot as plt
from PIL import Image

def neutron(N,ata,position):
    i=0
    while i<N:
        neutron=[]
        t=[0.0]
        v=7
        vitesse=[rd.uniform(-1*v,v),rd.uniform(-1*v,v),rd.uniform(-1*v,v)]
        position=position.copy()
        neutron.append(position)
        neutron.append(t)
        neutron.append(vitesse)
        ata.append(neutron)
        i=i+1
    return ata

def neutronvect(N, ata0, ata1, ata2, position):
    for i in range(N):
        t = [0.0]
        v = 7
        vitesse = [rd.uniform(-v, v), rd.uniform(-v, v), rd.uniform(-v, v)]
        ata0 = np.append(ata0, [position], axis=0)
        ata1 = np.append(ata1, [vitesse], axis=0)
        ata2 = np.append(ata2, [t], axis=0)
    return (ata0, ata1, ata2)

@njit
def desintegration_aux(ata0, ata2):
    dmin = 3
    for i in range(len(ata0)):
        dif = ata0 - ata0[i]
        distances = np.sqrt(np.sum(dif*dif, axis=1))
        nearby_indices = np.where(distances <= dmin)[0]
        bulle = 0
        for j in range(len(nearby_indices)):
            if i == j:
                continue
            bulle += 1
            if bulle > 50:
                ata2[nearby_indices[j]][0] = np.inf

def desintegration(ata0, ata1, ata2, dinf):
    age_fermi = 2
    desintegration_aux(ata0, ata2)
    indice = [i for i in range(len(ata0)) if ata2[i][0] <= age_fermi and np.sqrt(np.sum(ata0[i]**2)) < dinf]
    ata0 = np.array([ata0[i] for i in indice])
    ata1 = np.array([ata1[i] for i in indice])
    ata2 = np.array([ata2[i] for i in indice])
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

def fission(ata0, ata1, ata2, comb, r):
    i=0
    E_th=1.8
    v=7
    n=len(ata0)
    points = comb
    while i<n:
        a = ata0[i]
        distances = np.linalg.norm(points - a,axis=1)
        nearby_indices = np.where(distances <= r)[0]
        if nearby_indices.size == 0:
            i += 1
            continue
        noyau = np.argmin(distances)
        Ec=np.sqrt(np.sum(ata1[i]**2))
        P = (Ec - E_th) / (Ec * (1 + Ec / (2 * E_th)))
        if(rd.normal(0,1)<P):
            ata0, ata1, ata2 = neutronvect(3, ata0, ata1, ata2, points[noyau])
            ata2[i][0]=np.inf
        else:
            ata1[i]=-ata1[i]*10**(ata2[i][0])
            ata2[i][0]=0
        i=i+1
    return ata0, ata1, ata2

def combustible(r,d,n,e,rep):
    comb=[]
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
                    comb.append([xf,yf,j])
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
        image_files.append(os.path.join(image_folder,str(n)+'o.png'))
        n=n+1
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
    clip.write_videofile(image_folder+'ball-o.mp4')
    n=1
    while n<=i:
        os.remove('./'+dossier+"/"+str(n)+"o.png")
        n=n+1

def cycle(n,t,r,comb,ata,dinf,d,dossier):
    start = datetime.now()
    i=0
    tmp=len(ata)
    kmoy=0
    keff=[]
    neutron=[]
    neutron.append(len(ata))
    comb = np.array(comb)
    ata0 = np.array([x[0] for x in ata])
    ata1 = np.array([x[2] for x in ata])
    ata2 = np.array([x[1] for x in ata])
    while i<n:
        ata0,ata2=deplacement(ata0,ata1,ata2,t)
        ata0,ata1,ata2=fission(ata0,ata1,ata2,comb,r)
        ata1=moderation(ata1,t)
        ata0,ata1,ata2=desintegration(ata0, ata1, ata2, dinf)
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
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        draw_points1(comb)
        draw_points(ata0)
        pygame.display.flip()
        pygame.image.save(f, './'+dossier+"/"+str(i)+"o.png")
    delta=datetime.now() - start
    print ("Temps d'exécution : ", delta)
    movie(i,dossier)
    fig = plt.figure()
    plt.plot([i for i in range(len(neutron))],neutron)
    plt.title("Variation de Neutrons",fontsize = 14)
    plt.savefig('./'+dossier+'/var-o.png')
    plt.close('all')
    fig = plt.figure()
    plt.plot([i for i in range(len(keff))],keff)
    plt.title("Variation de keff",fontsize = 14)
    plt.savefig('./'+dossier+'/keff-o.png')
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
    cycle(n,t,uma,combustible(r,d,b,e,rep),neutron(100,[],[0,0,d/2]),r*rep,d,dossier)

pygame.init()
display = (800, 600)
f = pygame.display.set_mode(display, DOUBLEBUF | OPENGL)

# Paramètres de la vue 3D
gluPerspective(80, (display[0] / display[1]), 0.1, 150.0)
glTranslatef(0.0, 0.0, -40)

# Fonction pour dessiner les points
def draw_points(points):
    glBegin(GL_POINTS)
    for point in points:
        glColor3f(0.0, 0.0, 1.0)
        glVertex3fv(point)
    glEnd()

def draw_points1(points):
    glBegin(GL_POINTS)
    for point in points:
        glColor3f(0.0, 1.0, 0.0)
        glVertex3fv(point)
    glEnd()

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