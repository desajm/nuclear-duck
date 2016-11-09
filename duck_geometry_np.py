# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 19:41:33 2015

@author: sthagon
"""

from duck import spher_2_cart
from numpy import deg2rad, cross, dot, subtract, array
from numpy.linalg import norm

DEBUG = False


"""
ideja
OPIS DETEKTORA
- detektor definiramo s R, theta i phi i cinjenicom da je okomit na normalu
- onda to prebacimo u kartezijev sustav gdje je definiran s tockom i normalom
- tocku direktno dobijemo prebacivanjem R,theta,phi u x,y,z
- normalu dobijemo iz te tocke i ishodista jer je detektor okomit na ishodiste 
- x, y i z postavljeni su kao na  slici slika_osi.png
- theta ide od 0 do 180, a phi od 0 do 180 i 0 do "-180" s donje strane

OPIS CESTICE
- cesticu opisujemo pravcem
- krece iz ishodista i smjer joj je definiran s ishodistom i r,theta,phi -> x,y,z

SJECISTE
- sympy nam izbacuje sjeciste ravnine kojom je opisan detektor i pravca kojim opisujemo putanju cestice
- rezultat je tocka ako se sijeku
- sada je jos potrebno izracunati da li je to unutar detektora ili ne i to unutar kojeg pixela

UNUTAR DETEKTORA
- prebacimo se iz 3D sustava u 2D sustav ravnine
- http://stackoverflow.com/questions/23472048/projecting-3d-points-to-2d-plane
- imamo normalu i tocku ishodista
- treba definirati kooridnatni sustav 2D ravnine
- y' os treba biti ista kao i u 3D sustavu
- a x' os treba biti definirana na nacin da bude okomita na y' i normalu
"""

""" definiramo ravninu detektora
s normalom - usmjerenjem ravnine
i originom - tockom na ravnini
s obzirom na nasu situaciju mozemo ih definirati kao dva identicna vektora
"""
def create_detector_plane(R,theta,phi):    
    x,y,z = spher_2_cart(R,theta,phi)
    
    origin = array([x,y,z])    
    normal = array([x,y,z])
    
    return origin, normal


def create_particle_line(theta,phi):
    """ odredujemo dva vektora
    jedan je ishodiste (tocka na pravcu)
    drugi je usmjerenje sa theta i phi,
    a nalazi se na proizvoljnoj udaljenosti 1
    """
    point = array([0,0,0])
    direction = array(spher_2_cart(1.,theta,phi))
    
    return point, direction


def create_2D_system(normal):
    new_y_axis = [0,1,0]
    new_x_axis = cross(new_y_axis, normal)
    new_x_axis = new_x_axis / norm(new_x_axis)
    
    return new_x_axis, new_y_axis
    
def calculate_2D_coordinates(point, origin, x_axis, y_axis):
    helper = subtract(point,origin)
    #print helper, x_axis, y_axis
    return dot(x_axis,helper), dot(y_axis,helper)


def get_plane_impact_coordinates(intersect, detector):

    
    """ sjeciste u koordinatnom sustavu detektora
    ovisno o visini i sirini detektora
    sjeciste je unutar povrsine detektora ili ne
    """
    x, y = calculate_2D_coordinates(
                    intersect, 
                    detector['origin'], 
                    detector['cs_x'], 
                    detector['cs_y']
                )
    #print "impact coordinates on sensor"
    #print "x: ", x
    #print "y: ", y
    
    return x, y
    

def calculate_is_detected(detector, x, y):    
    width = detector['width']
    height = detector['height']   
    
    if (width/2.) > abs(x) and (height/2.) > abs(y):
        return True
    else:
        return False


def is_particle_geometrically_detected(detector, particle):
    """ sjeciste ravnine i pravca
    kad bi detektor imao beskonacnu povrsinu
    onda bi u ovoj tocki cestica pogadala detektor
    ako nema sjecista to znaci da cestica ide u paraleli s detektorm
    """
        
    ndotu = detector['normal'].dot(particle['direction'])
    
    if DEBUG:
        print detector['normal']
        print particle['direction']
        print ndotu
    
    """ ako je ndotu jednak nuli onda nema sjecista
    """
    epsilon=1e-6    
    if abs(ndotu) < epsilon:
        return False
    
    w = particle['point'] - detector['origin']
    si = -detector['normal'].dot(w) / ndotu
    
    intersect = w + si * particle['direction'] + detector['origin']
           
    if len(intersect) == 0:
        return False
 
    """ provjera da li cestica ide prema meti
    ili mozda od mete
    """
    p1 = intersect
    p2 = particle['direction']
    for i in range(3):
        if not p1[i] * p2[i] >= 0:
            return False
    
    #print "intersection: ", intersect, particle['line']
    
    """ dohvat koordinate na ravnini detektora
    potrebno je prvo postaviti koordinatni sustav unutar detektora
    """
    x, y = get_plane_impact_coordinates(intersect, detector)
    
    if DEBUG:
        print "X, Y: ", x, y

    """ izracun da li je unutra ili vani
    s obzirom na dimenzije detektora
    """        
    return calculate_is_detected(detector, x, y)
    

""" popunjavamo detektor polje
s potrebnim dodatnim podacima
"""
def initialize_detector(detector):
    """ izracun ravnine detektora
    racunamo ravninu u kojoj se nalazi detektor
    """
    detector['origin'], detector['normal'] = create_detector_plane(
                            detector['distance'],
                            detector['theta'],
                            detector['phi']
                            )                        
    
    """ koordinatni sustav detektora
    kreiramo 2D koordinatni sustav detektorove ravnine
    """
    detector['cs_x'], detector['cs_y'] = create_2D_system(detector['normal'])
    
    if DEBUG:
        print "2d system: ", detector['cs_x'], detector['cs_y']

    return detector    


""" popunjavamo particle dictionary
s formulom pravca
"""
def initialize_particle(particle):
    particle['point'], particle['direction'] = create_particle_line(
                        particle['theta'],
                        particle['phi']
                        )
    return particle



""" ovo koristimo za testiranje
pozovemo skriptu sa python duck_geometry.py, 
a prije toga popunimo ovdje odgovarajuce parametre detektora i cestice
"""

if __name__ == '__main__':
    
    """ podaci o detektoru
    theta_1 i phi_1 su kutevi
    distance je radius, tj. udaljenost od izvora
    """
    detector = {}
    detector['theta'] = deg2rad(50.0)
    detector['phi'] = deg2rad(0.)
    detector['distance'] = 0.05975
    detector['width'] = 0.05
    detector['height'] = 0.05
    
    initialize_detector(detector)    
    
    
    for i in range(100):
        """ podaci o cestici
        theta_2 i phi_2 definiraju pravac u kojem leti cestica
        """
        theta = 0.+ i
        phi = 0.
        particle = {}
        particle['theta'] = deg2rad(theta)
        particle['phi'] = deg2rad(phi)
        
        initialize_particle(particle)
    
        print theta, phi, 
        
        is_detected = is_particle_geometrically_detected(detector, particle)
    
        print is_detected












