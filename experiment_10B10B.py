DEPENDENCY_MODULE = 'duck_mc'

E_12_START = 0. 
E_12_STOP = 30.0
E_12_STEP = 0.1

BROJ_PONAVLJANJA = 10

def duck_mc_experiment(E_12, broj_ponavljanja, postav_index):
    """ postavke detektora na BoBo eksperimentu
    """
    POSTAV = {
        'postav_1': {'distance': [0.35875,0.35675,0.35775,0.35975], 'theta': [40.49, 20.00, 19.995, 40.095], 'phi': [0, 0, 180, 180], 'num_stripes': 16, 'detector_width': 0.05},
        'postav_2': {'distance': [0.35875,0.35675,0.35775,0.35975], 'theta': [40.49, 20.00, 30.00, 50.00], 'phi': [0, 0, 180, 180], 'num_stripes': 16, 'detector_width': 0.05},
        'postav_3': {'distance': [0.35875,0.35675,0.35775,0.35975], 'theta': [46.49, 26.00, 33.00, 53.00], 'phi': [0, 0, 180, 180], 'num_stripes': 16, 'detector_width': 0.05}
        }    
        

    kwargs = {
        "E_p": 72.132,
        "m_p": 10.,
        "m_t": 10.,
        "m_X": 12.,
        "m_E": 8.,
        "m_C": 10.,
        "m_D": 2.,
        "Q_01": 19.15971, # Q_0 za reakciju 1
        "E_prag": 25.1868,
        "Emin_C": [37.4, 38.6, 35.4, 34.5], #za 10B za det na 0st
        "Emin_D": [3.2, 3.3, 3.0, 3.0] #za d, za det na 0 st
        
    }


    detectors = [duck_mc.define_detector(POSTAV, i, postav_index) for i in range(4)]
        
    try:
        result = duck_mc.calculate_event_detection_E12(detectors, E_12, broj_ponavljanja, **kwargs)
        #result = duck_14N_test.calculate_event_detection_E12(detectors, E_12, BROJ_PONAVLJANJA, E_prag, Q_01)
    except:
        e = sys.exc_info()[0]
        return e
    
    return result


