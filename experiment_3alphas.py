DEPENDENCY_MODULE = 'duck_mc_3alphas'

E_12_START = 0. 
E_12_STOP = 30.0
E_12_STEP = 1.

BROJ_PONAVLJANJA = 1000

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


    detectors = [duck_mc_3alphas.define_detector(POSTAV, i, postav_index) for i in range(4)]
        
    try:
        result = duck_mc_3alphas.calculate_event_detection_E12(detectors, E_12, broj_ponavljanja, **kwargs)
    except:
        e = sys.exc_info()[0]
        return e
    
    return result


