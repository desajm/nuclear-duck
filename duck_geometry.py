# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 19:41:33 2015

@author: sthagon
"""


from sympy import Point3D, Plane, Line3D, N
from duck import spher_2_cart
from numpy import deg2rad, cross, dot, subtract
from numpy.linalg import norm


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

def create_detector_plane(R,theta,phi):    
    x,y,z = spher_2_cart(R,theta,phi)
    
    origin = Point3D(x,y,z)
    
    normal_vector = (x,y,z)
    
    return Plane(origin, normal_vector), origin, normal_vector


def create_particle_line(theta,phi):
    """ odredujemo dvije tocke
    jedna je ishodiste
    druga je usmjerena sa theta i phi,
    a nalazi se na proizvoljnoj udaljenosti 1
    """
    p_0 = Point3D(0,0,0)
    p_1 = Point3D(spher_2_cart(1.,theta,phi))
    
    return Line3D(p_0,p_1)


def create_2D_system(normal):
    new_y_axis = [0,1,0]
    new_x_axis = cross(new_y_axis, normal)
    new_x_axis = new_x_axis / norm(new_x_axis)
    #print "y_axis: ", new_y_axis
    #print "x_axis: ", new_x_axis
    
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
                    intersect[0], 
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
    intersect = detector['plane'].intersection(particle['line'])

    if len(intersect) == 0:
        return False
 
    """ provjera da li cestica ide prema meti
    ili mozda od mete
    """
    p1 = intersect[0].evalf()
    p2 = particle['line'].p2.evalf()
    for i in range(3):
        if not p1[i] * p2[i] >= 0:
            return False
    
    #print "intersection: ", intersect, particle['line']
    
    """ dohvat koordinate na ravnini detektora
    potrebno je prvo postaviti koordinatni sustav unutar detektora
    """
    x, y = get_plane_impact_coordinates(intersect, detector)
    
#    print x, y,

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
    detector['plane'], detector['origin'], detector['normal'] = create_detector_plane(
                            detector['distance'],
                            detector['theta'],
                            detector['phi']
                            )                        
    
    """ koordinatni sustav detektora
    kreiramo 2D koordinatni sustav detektorove ravnine
    """
    detector['cs_x'], detector['cs_y'] = create_2D_system(detector['normal'])

    return detector    


""" popunjavamo particle dictionary
s formulom pravca
"""
def initialize_particle(particle):
    particle['line'] = create_particle_line(
                        particle['theta'],
                        particle['phi']
                        )
    return particle


def calculate_intersection(plane, line):
    return plane.intersection(line)

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped



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











