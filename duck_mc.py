# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:35:29 2016

@author: desa
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 21 21:27:25 2015

@author: desa
"""
#from ROOT import *
import numpy as np
from duck import spher_2_cart, cart_2_spher
from duck_geometry_np import initialize_detector, initialize_particle, is_particle_geometrically_detected

"""
hist1 = TH2F('hist1', 'E(10B) - theta(10B)' ,1024, -10, 190, 1024, -10, 90)
hist1.SetOption("COLZ")
hist2 = TH2F('hist2', 'E(6Li) - theta(6Li)' ,1024, -10, 190, 1024, -10, 90)
hist2.SetOption("COLZ")
hist3 = TH2F('hist3', 'phi(6Li) - phi(4He)', 1024, -200, 200, 1024, -200, 200)
hist3.SetOption("COLZ")
hist4 = TH2F('hist3', 'theta(10B) - theta(10B)', 1024, -20, 200, 1024, -20, 200)
hist4.SetOption("COLZ")
hist5 = TH2F('hist1', 'E(6Li) - E(4He)' ,1024, -10, 90, 1024, -10, 90)
hist5.SetOption("COLZ")
"""



def define_detector(POSTAV, i, postav_key):
    detector = {}
    detector['distance'] = POSTAV[postav_key]['distance'][i]
    detector['theta'] = np.deg2rad(POSTAV[postav_key]['theta'][i])
    detector['phi'] = np.deg2rad(POSTAV[postav_key]['phi'][i])
    detector['num_stripes'] = POSTAV[postav_key]['num_stripes']
    detector['width'] = POSTAV[postav_key]['detector_width']
    detector['height'] = detector['width']

    initialize_detector(detector)
    
    return detector


# phi ide od -180 do 180
def get_random_phi():
    return np.deg2rad(360. * np.random.random_sample() - 180.)

# theta ide od 0. do 180.
def get_random_theta():
    cos_theta = 2. * np.random.random_sample() - 1.
    return np.arccos(cos_theta)
    
def is_condition_20(theta, m_1, m_2, m_p, Q, E_p):
    cos_2_theta = np.cos(theta)**2
    temp_1 = m_2 * (m_2 + m_1) / (m_p * m_1)
    temp_2 = m_p / m_2 - 1 - Q / E_p
    if cos_2_theta >= temp_1 * temp_2:
        return True
    return False


""" Racuna da li je meta mogla detektirati cesticu
dva uvjeta:
- cestica je morala ici u smjeru mete
- cestica je morala imati dovoljno energije
output:
- ind_detector ako je cestica detektirana
- -1 ako nije detektirana
"""
def is_particle_detected(particle, detectors, min_energies):
    for ind_detector in range(len(detectors)):
        detector = detectors[ind_detector]
        if is_particle_geometrically_detected(detector, particle):
            if particle['E_lab'] >= min_energies[ind_detector]:
                return ind_detector
            else:
                return -1
    return -1    


def is_2particle_event_detected(detectors, theta_X, phi_X, Q_1, E_12, E_p, m_p, m_t, m_X, m_E, m_C, m_D, Emin_C, Emin_D, **kwargs):
    
    E_in = E_p * m_t / (m_t + m_p)
    
    #print "E_in", E_in
        
    E_out = E_in  + Q_1
        
    #print "E_out", E_out
    
    E_X = E_out * m_E / (m_X + m_E)
    
        
    E_E = E_out - E_X
    
    
    theta_E = np.deg2rad(180 - np.rad2deg(theta_X))
            
    if np.rad2deg(phi_X) > 0.0:
        phi_E = np.deg2rad(np.rad2deg(phi_X)-180)
            
    if np.rad2deg(phi_X) < 0.0:
        phi_E = np.deg2rad(np.rad2deg(phi_X) + 180)
     
     
    #print "E_X, E_E", E_X, E_E
            
    """
    treba nam iznos brzine čestice X u laboratorijskom sustavu. Računamo ga iz E_X_lab
    """   
        
    a_X = np.sqrt(m_X * m_p * E_p)/(m_t + m_p)
    E_X_lab = E_X + 2.0 * np.sqrt(E_X) * a_X * np.cos(theta_X) + a_X * a_X
    v_X_iznos_lab = np.sqrt(2.0 * E_X_lab / m_X)
        
    theta_X_lab =  np.arccos((np.sqrt(E_X) * np.cos(theta_X) + a_X) / (np.sqrt(E_X + 2.0 * np.sqrt(E_X) * a_X * np.cos(theta_X) + a_X * a_X)))       
    phi_X_lab = phi_X
        
    v_X_cart_lab = spher_2_cart(v_X_iznos_lab, theta_X_lab, phi_X_lab)
      
    #print "thetaX_lab",  np.rad2deg(theta_X_lab)  
    #print  v_X_cart_lab   
       
      
      
    a_E = np.sqrt(m_E * m_p * E_p)/(m_t + m_p)
    E_E_lab = E_E + 2.0 * np.sqrt(E_E) * a_E * np.cos(theta_E) + a_E * a_E
    v_E_iznos_lab = np.sqrt(2.0 * E_E_lab / m_E)
        
    theta_E_lab =  np.arccos((np.sqrt(E_E) * np.cos(theta_E) + a_E) / (np.sqrt(E_E + 2.0 * np.sqrt(E_E) * a_E * np.cos(theta_E) + a_E * a_E)))       
    phi_E_lab = phi_E
    
    #hist4.Fill(np.rad2deg(theta_E_lab), np.rad2deg(theta_X_lab)) 
    """
    RASPAD X -> C + D
    slucajno odabiremo theta, phi za cesticu C
    kuteve biramo u odnosu na smjer gibanja cestice X
    """ 
    theta_C = get_random_theta()
    phi_C = get_random_phi()
        
    #print  "theta_C, phi_C", np.rad2deg(theta_C), np.rad2deg(phi_C)
    """
    Energija koja je na raspolaganju za ovaj raspad je E_12
    """        
    E_C = m_D * E_12 / (m_C + m_D)
        
        
    """ brzina cestice C
    """
    v_C_iznos = np.sqrt(2.0 * E_C / m_C)
    v_C_cart = spher_2_cart(v_C_iznos, theta_C, phi_C)
    
    #print "v_C_cart", v_C_cart
        
        
    """ brzina čestice D
    """
    E_D = m_C * E_12 / (m_C + m_D)
        
    theta_D = np.deg2rad(180 - np.rad2deg(theta_C))
            
    if np.rad2deg(phi_C) > 0.0:
        phi_D = np.deg2rad(np.rad2deg(phi_C)-180)
            
    if np.rad2deg(phi_C) < 0.0:
        phi_D = np.deg2rad(np.rad2deg(phi_C) + 180)
        
    v_D_iznos = np.sqrt(2 * E_D / m_D)
    v_D_cart = spher_2_cart(v_D_iznos, theta_D, phi_D)
    
    #print "v_D_cart", v_D_cart
    
    #print  "theta_D, phi_D", np.rad2deg(theta_D), np.rad2deg(phi_D)
    """
    TRANSFORMACIJA IZ SUSTAVA U KOJEM X MIRUJE U LABORATORIJSKI SUSTAV
    """
    v_C_lab_cart = v_C_cart[0] + v_X_cart_lab[0], v_C_cart[1] + v_X_cart_lab[1], v_C_cart[2] + v_X_cart_lab[2]
                
    v_C_lab_iznos = np.sqrt(v_C_lab_cart[0] * v_C_lab_cart[0] + v_C_lab_cart[1] * v_C_lab_cart[1] + v_C_lab_cart[2] * v_C_lab_cart[2])
            
    E_C_lab = m_C * v_C_lab_iznos * v_C_lab_iznos/2.0  
        
    v_C_lab_spher =  cart_2_spher(v_C_lab_cart[0], v_C_lab_cart[1], v_C_lab_cart[2])      
        
    #print "thetaC_lab, phiC_lab", np.rad2deg(v_C_lab_spher[1]), np.rad2deg(v_C_lab_spher[2])
    #hist2.Fill(np.rad2deg(v_C_lab_spher[1]), E_C_lab)  
    
    """
    PARAMETRI KOJI NAM TREBAJU SU E_C_lab, theta_C_lab, phi_C_lab
                                      E_D_lab, theta_D_lab, phi_D_lab    
    theta_C_lab = np.rad2deg(v_C_lab_spher[1])
    phi_C_lab = np.rad2deg(v_C_lab_spher[2]) 
    """
      
    v_D_lab_cart = v_D_cart[0] + v_X_cart_lab[0], v_D_cart[1] + v_X_cart_lab[1], v_D_cart[2] + v_X_cart_lab[2]
        
        
    v_D_lab_iznos = np.sqrt(v_D_lab_cart[0] * v_D_lab_cart[0] + v_D_lab_cart[1] * v_D_lab_cart[1] + v_D_lab_cart[2] * v_D_lab_cart[2])
            
    E_D_lab = m_D * v_D_lab_iznos * v_D_lab_iznos/2.0   
            
    v_D_lab_spher =  cart_2_spher(v_D_lab_cart[0], v_D_lab_cart[1], v_D_lab_cart[2]) 

   
    #print "thetaD_lab, phiD_lab", np.rad2deg(v_D_lab_spher[1]), np.rad2deg(v_D_lab_spher[2])     
    
    #hist1.Fill(np.rad2deg(theta_X_lab), E_X_lab)
    #hist3.Fill(np.rad2deg(v_D_lab_spher[2]), np.rad2deg(v_C_lab_spher[2])) 
    #hist3.Fill(np.rad2deg(phi_X_lab), np.rad2deg(phi_E_lab))
    
    #hist5.Fill(E_D_lab, E_C_lab)  
 
    """
    Prvo gledamo je li detektirana čestica C jer je to 6Li, a njih ima manje nego 4He
    Ako jest, onda gledamo je li detektirana 4He
    """
        
    particle_C = {}
    particle_C['theta'] = v_C_lab_spher[1]
    particle_C['phi'] = v_C_lab_spher[2]
    particle_C['E_lab'] = E_C_lab
    initialize_particle(particle_C)
        
    ind_C_detected = is_particle_detected(particle_C, detectors, Emin_C)
        
    """ ako cestica C nije detektirana
    onda odmah skoci na slijedeci event
    """
    if ind_C_detected < 0:
        return False
        
        
    particle_D = {}
    particle_D['theta'] = v_D_lab_spher[1]
    particle_D['phi'] = v_D_lab_spher[2]
    particle_D['E_lab'] = E_D_lab
    initialize_particle(particle_D)
                
                
    ind_D_detected = is_particle_detected(particle_D, detectors, Emin_D)
        
    """ ako cestica D nije detektirana
    onda odmah skoci na slijedeci event
    """
    if ind_D_detected < 0:
        return False            
    
    
    
   
    result = [ind_C_detected, ind_D_detected]
    
    return result


def calculate_event_detection_E12(detectors, E_12, BROJ_PONAVLJANJA, **kwargs):
    """ resetiramo brojac na 0 """
    num_det = 4
    broj_detektiranih_cestica = [[0 for x in range(num_det)] for y in range(num_det)] 

    """
    ZA SVAKI E_12 IZVRTIMO PAR TISUĆA DOGAĐAJA (slučajnih smjerova X)
    """
    brojac = 0

    while brojac < BROJ_PONAVLJANJA:
        brojac += 1
        
        """
        IZRACUNAVAMO EFEKTIVNU Q-VRIJEDNOST S KOJOM RACUNAMO KINEMATIKU
        E_prag - energija praga jezgre X za raspad u kanal C + D
        Q0 - zadajemo ovisno o reakciji koju zelimo promatrati
        """
        E_pob = kwargs['E_prag'] + E_12
        Q_1 = kwargs['Q_01'] - E_pob

        
        """
        dohvacamo random thetu i phi za cesticu X, u radijanima
        kod thete generiramo random cos_theta uniformnom distribucijom
        phi generiramo direktno

        """
            
        theta_X = get_random_theta()
        phi_X = get_random_phi()
        
        #print ""
        #print "theta_X, phi_X", np.rad2deg(theta_X), np.rad2deg(phi_X)
        

        """ ako su i C i D detektirane
        onda povecamo brojac detektiranih cestica
        """
        ind_detected = is_2particle_event_detected(detectors, theta_X, phi_X, Q_1, E_12, **kwargs)
        #ind_detected = is_2particle_event_detected(detectors, theta_X, phi_X, Q_1, E_12, E_p, m_p, m_t, m_X, m_E, m_C, m_D, Emin_C, Emin_D)
        if not ind_detected:
            continue

        """ povecavamo statistiku brojaca
        za neku od kombinacija detektora """
        broj_detektiranih_cestica[ind_detected[0]][ind_detected[1]] += 1


    """ logiramo broj detektiranih cestica za zadani E_pob
    """
    log_line = "statistika, %f, %f" % (E_12, E_pob)    
    for i in range(num_det):
        for j in range(num_det):
            log_line += ", %d"%(broj_detektiranih_cestica[i][j])
    print log_line
    return log_line





"""
Promatramo reakciju p + t -> X + E (gdje E ostaje nedetektirana) a X se raspada dalje na X - > C + D
Sve skupa imamo p + t -> C + D + E gdje E nedetektiramo, a C i D detektiramo
E_12 - energija iznad praga koji je potreban za raspad X -> C + D
E_MAX_REL - maksimalna vrijednost do koje E12 može ići, odgovara maksimalnoj 
relativnoj energiji koju vidimo u spektru za jezgru X.
"""

if __name__ == '__main__':

    POSTAV = {
        'postav_1': {'distance': [0.35875,0.35675,0.35775,0.35975], 'theta': [40.49, 20.00, 19.995, 40.095], 'phi': [0, 0, 180, 180], 'num_stripes': 16, 'detector_width': 0.05},
        'postav_2': {'distance': [0.35875,0.35675,0.35775,0.35975], 'theta': [40.49, 20.00, 30.00, 50.00], 'phi': [0, 0, 180, 180], 'num_stripes': 16, 'detector_width': 0.05},
        'postav_3': {'distance': [0.35875,0.35675,0.35775,0.35975], 'theta': [46.49, 26.00, 33.00, 53.00], 'phi': [0, 0, 180, 180], 'num_stripes': 16, 'detector_width': 0.05}
        }
    
        
    """
    PRVI SLUČAJ: 10B + 10B -> 10B + 10B -> 6Li + 4He + 10B
    """
    
    
    
    
    postav_index = 'postav_1'
    detectors = [define_detector(POSTAV, i, postav_index) for i in range(4)]
 
    kwargs = {
        "E_p": 72.132,
        "m_p": 10.,
        "m_t": 10.,
        "m_X": 10.,
        "m_E": 10.,
        "m_C": 6.,
        "m_D": 4.,
        "Q_01": 0., # Q_0 za reakciju 1
        "E_prag": 4.461,
        "Emin_C": [17.5, 18.0, 16.6, 16.2], #za det na 0st
        "Emin_D": [9.4, 9.7, 8.9, 8.7],
        #"Emin_C": [21.6, 19.0, 17.5, 19.9], #za det na najvisem kutu postav 1 (rub detektora)
        #"Emin_D": [11.5, 10.1, 9.3, 10.6], #za det na najvisem kutu postav 1 (rub detektora)
       # "Emin_C": [20.0, 18.3, 16.8, 18.5], #za det na najmanjem kutu postav 1 (rub detektora)
        #"Emin_D": [10.7, 9.8, 9.0, 9.9], #za det na najmanjem kutu postav 1 (rub detektora)
    }

   
    
    BROJ_PONAVLJANJA = 10000
    E_MAX_REL = 30.0
    E_12_STEP = 10
    E_12 = 0. 
    


    while E_12 <= E_MAX_REL:
        calculate_event_detection_E12(detectors, E_12, BROJ_PONAVLJANJA, **kwargs)
        E_12 += E_12_STEP
    
    
