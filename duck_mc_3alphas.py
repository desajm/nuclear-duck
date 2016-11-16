# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime
from duck import spher_2_cart, cart_2_spher
from duck_geometry_np import initialize_detector, initialize_particle, is_particle_geometrically_detected

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


def is_3particle_event_detected(detectors, theta_X, phi_X, Q_1, E_12, E_p, m_p, m_t, m_X, m_E, m_C, m_D, m_F, m_G, Emin_F, Emin_G, Emin_D, **kwargs):
    
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
    treba nam iznos brzine čestice X u laboratorijskom sustavu. Racunamo ga iz E_X_lab
    """   
        
    a_X = np.sqrt(m_X * m_p * E_p)/(m_t + m_p)
    E_X_lab = E_X + 2.0 * np.sqrt(E_X) * a_X * np.cos(theta_X) + a_X * a_X
    v_X_iznos_lab = np.sqrt(2.0 * E_X_lab / m_X)
        
    theta_X_lab =  np.arccos((np.sqrt(E_X) * np.cos(theta_X) + a_X) / (np.sqrt(E_X + 2.0 * np.sqrt(E_X) * a_X * np.cos(theta_X) + a_X * a_X)))       
    phi_X_lab = phi_X
        
    v_X_cart_lab = spher_2_cart(v_X_iznos_lab, theta_X_lab, phi_X_lab)
      
      
      
    a_E = np.sqrt(m_E * m_p * E_p)/(m_t + m_p)
    E_E_lab = E_E + 2.0 * np.sqrt(E_E) * a_E * np.cos(theta_E) + a_E * a_E
    v_E_iznos_lab = np.sqrt(2.0 * E_E_lab / m_E)
        
    theta_E_lab =  np.arccos((np.sqrt(E_E) * np.cos(theta_E) + a_E) / (np.sqrt(E_E + 2.0 * np.sqrt(E_E) * a_E * np.cos(theta_E) + a_E * a_E)))       
    phi_E_lab = phi_E
    
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
        
    
    
      
    v_D_lab_cart = v_D_cart[0] + v_X_cart_lab[0], v_D_cart[1] + v_X_cart_lab[1], v_D_cart[2] + v_X_cart_lab[2]
        
        
    v_D_lab_iznos = np.sqrt(v_D_lab_cart[0] * v_D_lab_cart[0] + v_D_lab_cart[1] * v_D_lab_cart[1] + v_D_lab_cart[2] * v_D_lab_cart[2])
            
    E_D_lab = m_D * v_D_lab_iznos * v_D_lab_iznos/2.0   
            
    v_D_lab_spher =  cart_2_spher(v_D_lab_cart[0], v_D_lab_cart[1], v_D_lab_cart[2]) 

    """
    Prvo gledamo je li detektirana čestica D, ako jest idemo na sljedeci raspad
    """
        
    particle_D = {}
    particle_D['theta'] = v_D_lab_spher[1]
    particle_D['phi'] = v_D_lab_spher[2]
    particle_D['E_lab'] = E_D_lab
    initialize_particle(particle_D)
        
    ind_D_detected = is_particle_detected(particle_D, detectors, Emin_D)
        
    """ ako cestica C nije detektirana
    onda odmah skoci na slijedeci event
    """
    if ind_D_detected < 0:
        return False
    
    """
    RASPAD C -> F + G
    slucajno odabiremo theta, phi za cesticu F
    kuteve biramo u odnosu na smjer gibanja cestice C
    """ 
    theta_F = get_random_theta()
    phi_F = get_random_phi()
    
    #print  theta_C, phi_C
    """
    Energija koja je na raspolaganju za ovaj raspad je UP

    """        
    
    E_UP = 0.09184
    E_F = m_G * E_UP / (m_F + m_G)
    
    
    """ brzina cestice F
    """
    v_F_iznos = np.sqrt(2.0 * E_F / m_F)
    v_F_cart = spher_2_cart(v_F_iznos, theta_F, phi_F)
    
    
    """ brzina čestice G
    """
    E_G = m_F * E_UP / (m_F + m_G)
    
    theta_G = np.deg2rad(180 - np.rad2deg(theta_F))
        
    if np.rad2deg(phi_F) > 0.0:
        phi_G = np.deg2rad(np.rad2deg(phi_F)-180)
        
    if np.rad2deg(phi_F) < 0.0:
        phi_G = np.deg2rad(np.rad2deg(phi_F) + 180)
    
    v_G_iznos = np.sqrt(2 * E_G / m_G)
    v_G_cart = spher_2_cart(v_G_iznos, theta_G, phi_G)


    """
    TRANSFORMACIJA IZ SUSTAVA U KOJEM C MIRUJE U LABORATORIJSKI SUSTAV
    """
    v_F_lab_cart = v_F_cart[0] + v_C_lab_cart[0], v_F_cart[1] + v_C_lab_cart[1], v_F_cart[2] + v_C_lab_cart[2]
            
    v_F_lab_iznos = np.sqrt(v_F_lab_cart[0] * v_F_lab_cart[0] + v_F_lab_cart[1] * v_F_lab_cart[1] + v_F_lab_cart[2] * v_F_lab_cart[2])
        
    E_F_lab = m_F * v_F_lab_iznos * v_F_lab_iznos/2.0  
    
    v_F_lab_spher =  cart_2_spher(v_F_lab_cart[0], v_F_lab_cart[1], v_F_lab_cart[2])   
    
    #print v_C_lab_iznos, v_C_lab_spher[0]
    
  
    v_G_lab_cart = v_G_cart[0] + v_C_lab_cart[0], v_G_cart[1] + v_C_lab_cart[1], v_G_cart[2] + v_C_lab_cart[2]
    
    v_G_lab_iznos = np.sqrt(v_G_lab_cart[0] * v_G_lab_cart[0] + v_G_lab_cart[1] * v_G_lab_cart[1] + v_G_lab_cart[2] * v_G_lab_cart[2])
        
    E_G_lab = m_G * v_G_lab_iznos * v_G_lab_iznos/2.0  
    
    v_G_lab_spher =  cart_2_spher(v_G_lab_cart[0], v_G_lab_cart[1], v_G_lab_cart[2]) 
    
    """
    GLEDAMO JE LI ČESTICE F I G POGAĐAJU DETEKTORE SLOŽENE U POSTAV 1 ILI 2 ILI 3
    """
    
    
    particle_F = {}
    particle_F['theta'] = v_F_lab_spher[1]
    particle_F['phi'] = v_F_lab_spher[2]
    particle_F['E_lab'] = E_F_lab
    initialize_particle(particle_F)
    
    ind_F_detected = is_particle_detected(particle_F, detectors, Emin_F)
    
    """ ako cestica F nije detektirana
    onda odmah skoci na slijedeci event
    """
    if ind_F_detected < 0:
        return False 
    
    
    particle_G = {}
    particle_G['theta'] = v_G_lab_spher[1]
    particle_G['phi'] = v_G_lab_spher[2]
    particle_G['E_lab'] = E_G_lab
    initialize_particle(particle_G)
            
            
    ind_G_detected = is_particle_detected(particle_G, detectors, Emin_G)
    
    """ ako cestica D nije detektirana
    onda odmah skoci na slijedeci event
    """
    if ind_G_detected < 0:
        return False 
    
   
    

    """ ispis u datoteku
    E_12, E_pob, brojac, theta_C, phi_C, theta_D, phi_D
    
    log_line = "%f, %f, %d, %f, %f, %f, %f" % (E_12,
                                               E_pob,
                                               brojac,
                                               np.rad2deg(particle_C['theta']),
                                               np.rad2deg(particle_C['phi']),
                                               np.rad2deg(particle_D['theta']),
                                               np.rad2deg(particle_D['phi']),
                                            )
    #print log_line
    target.write(log_line)
    target.write("\n")
    """
    result = [ind_D_detected, ind_F_detected, ind_G_detected]
    
    return result







    
    


def calculate_event_detection_E12(detectors, E_12, BROJ_PONAVLJANJA, **kwargs):
    """ resetiramo brojac na 0 """
    num_det = 4
    broj_detektiranih_cestica = [[[0 for x in range(num_det)] for y in range(num_det)] for z in range(num_det)]

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
        #ind_detected = is_3particle_event_detected(detectors, theta_X, phi_X, Q_1, E_12, E_p, m_p, m_t, m_X, m_E, m_C, m_D, m_F, m_G, Emin_F, Emin_G, Emin_D)
        ind_detected = is_3particle_event_detected(detectors, theta_X, phi_X, Q_1, E_12, **kwargs)
        if not ind_detected:
            continue

        """ povecavamo statistiku brojaca
        za neku od kombinacija detektora """
        broj_detektiranih_cestica[ind_detected[0]][ind_detected[1]][ind_detected[2]] += 1


    """ logiramo broj detektiranih cestica za zadani E_pob
    """
    log_line = "%s, %f, %f" % (str(datetime.now()), E_12, E_pob)    
    for i in range(num_det):
        for j in range(num_det):
            for k in range(num_det):
                log_line += ", %d"%(broj_detektiranih_cestica[i][j][k])
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
    
    postav_index = 'postav_1'
    detectors = [define_detector(POSTAV, i, postav_index) for i in range(4)]
    
    kwargs = {
            "E_p": 72.132,
            "m_p": 10.,
            "m_t": 10.,
            "m_X": 12.,
            "m_E": 8.,
            "m_C": 8.,
            "m_D": 4.,
            "m_F": 4.,
            "m_G": 4.,
            "Q_01": 16.12971, #(za vrh Q1, 19.15971 - 3.03) # Q_0 za reakciju 1
            "E_prag": 7.3666,
            "Emin_F": [9.4, 9.7, 8.9, 8.7], #za det na 0st
            "Emin_G": [9.4, 9.7, 8.9, 8.7],
            "Emin_D": [9.4, 9.7, 8.9, 8.7],
            }
    
    BROJ_PONAVLJANJA = 10000
    E_MAX_REL = 30.0
    E_12_STEP = 1.0
    E_12 = 0. 


    while E_12 <= E_MAX_REL:
        result = calculate_event_detection_E12(detectors, E_12, BROJ_PONAVLJANJA, **kwargs)
        print result
        E_12 += E_12_STEP
    
    
