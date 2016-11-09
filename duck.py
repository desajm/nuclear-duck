# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 21:56:04 2014

@author: desa
"""

import re, pexpect
from ROOT import TTree, TFile
from time import time
from datetime import timedelta
from joblib import Memory
import numpy as np
import math as m
from numpy.linalg import norm


mem = Memory(cachedir='cache', verbose=0)


def start_timer():
    """
    Return start time of the experiment.

    Example:
    start_time = start_timer()
    print_time(start_time, step, num_of_steps)
    """
    return time()


def print_time(start_time, step, num_of_steps, milestone=10000):
    """Print how much time elapsed and how much time is remaining till the end of experiment.
    
    Arguments:
    start_time -- time when experiment started (get with start_timer())
    step -- current step of the experiment
    num_of_steps -- total number of steps that need to be performed
    milestone -- after completing this number of steps print status
    
    Example:
    start_time = start_timer()
    for i in xrange(100000):
        print_time(start_time, i, 100000, 1000)
        do_something()
    """

    if step % milestone != 0:
        return
        
    elapsed_time = round(time() - start_time)
    percentage = (step * 1.) / num_of_steps
    remaining_time = 0
    
    if percentage != 0:
        remaining_time = round(elapsed_time * (1 - percentage) / percentage)
        remaining_time = timedelta(0, remaining_time)

    elapsed_time = timedelta(0, elapsed_time)
    print "%.1f%%, elapsed %s, remaining %s, event %d from %d " % (percentage*100, elapsed_time, remaining_time, step, num_of_steps)


def theil_sen(arr_x, arr_y):
    """
    Theil Sen is a method for robust linear regression that chooses the median slope among all lines through pairs of two-dimensional sample points.
    It's especially good when we don't have normally distributed noise.
    Problem is that it's pretty slow for big datasets.

    Arguments:
        arr_x -- array with x coordinates
        arr_y -- array with y coordinates
    
    Returns:
        m -- line slope (y = m*x +b)
        b -- line intercept of the line
        std_error -- standard error of regression
        r_squared -- coefficient of determination
    
    Example:
        import numpy as np
        x = np.random.rand(10)
        y = np.random.rand(10)
        slope,intercept,std_error,r_squared = theil_sen(x,y)
    """
    
    N = arr_x.shape[0]
    
    def slope(i, j):
        if arr_x[i] == arr_x[j]:
            return
        return float(arr_y[i] - arr_y[j]) / (arr_x[i] - arr_x[j])
 
    def median(L):
        L.sort()
        if len(L) & 1:
            return L[len(L) // 2]
        else:
            return (L[len(L) // 2 - 1] + L[len(L) // 2]) / 2.0

    tmp = [slope(i, j) for i in range(N) for j in range(i)]        
    med = median(tmp)
     
    def error(i):
        return arr_y[i] - med * arr_x[i]
 
    b = median([error(i) for i in range(N)])
    
    
    def calculate_r_squared():
        y_mean = np.mean(arr_y)
        ss_tot = sum((arr_y - y_mean)**2)
        f_back = med * arr_x + b
        ss_reg = sum((f_back - y_mean)**2)
        return 1 - ss_reg / ss_tot
    
    def calculate_std_error():
        f_y = med * arr_x + b
        N = arr_x.shape[0]
        return np.sqrt(sum((f_y - arr_y)**2) / N)
        
    std_error = calculate_std_error()
    r_squared = calculate_r_squared()
    
    return med, b, std_error, r_squared

def theil_sen_categories(data):
    num_categories = len(data)
    
    def slope(el_1, el_2):
        if el_1[0] == el_2[0]:
            return
        return float(el_1[1] - el_2[1]) / (el_1[0] - el_2[0])
 
    def median(L):
        L.sort()
        if len(L) & 1:
            return L[len(L) // 2]
        else:
            return (L[len(L) // 2 - 1] + L[len(L) // 2]) / 2.0
            
    tmp = []
    for i in xrange(num_categories):
        for j in xrange(i):
            for el_1 in data[i]:
                for el_2 in data[j]:
                    tmp.append(slope(el_1,el_2))
                   
    med = median(tmp)
     
    def error(el):
        return el[1] - med * el[0]
 
    tmp_err = []
    for category in data:
        for el in category:
            tmp_err.append(error(el))
        
    b = median(tmp_err)    
    
    return med, b


def stochastic_theil_sen_categories(data, sample_size=200, num_samples=20, uniform=False):
    
    arr_slope = np.zeros(num_samples)
    arr_intercept = np.zeros(num_samples)
        
    
    total_points = sum([len(item) for item in data])   

    num_sample_points = sample_size/len(data)
    


    for i in xrange(num_samples):
        sample = []
        for dataset in data:

            if len(dataset) < 10:
                continue

            if not uniform:
                num_sample_points = int(sample_size * (float(len(dataset))/total_points))

            #print 'category', num_sample_points, len(dataset)
                        
            population = np.random.randint(0, len(dataset), size = num_sample_points)
            
            tmp_points = []
            for j in xrange(num_sample_points):
                tmp_points.append(dataset[population[j]])

            sample.append(tmp_points)

        med,b = theil_sen_categories(sample)

        print i, med, b

        arr_slope[i] = med
        arr_intercept[i] = b
    
    slope = np.median(arr_slope)
    intercept = np.median(arr_intercept)
    
    return slope, intercept
        

def stochastic_theil_sen(arr_x, arr_y, sample_size=200, num_samples=20):
    """
    Calculate linear coefficients using Theil Sen method on several (num_samples) samples of sample_size several times. Return median of linear coefficients.
    
    This should be used when we have too many data points so that regular Theil Sen is too slow.    
    
    Theil Sen is a method for robust linear regression that chooses the median slope among all lines through pairs of two-dimensional sample points.
    It's especially good when we don't have normally distributed noise.
    Problem is that it's pretty slow for big datasets.

    Arguments:
        arr_x -- array with x coordinates
        arr_y -- array with y coordinates
        sample_size -- how many points in a sample (default = 200)
        num_samples -- how many samples we should use (default = 20)
        arr_type -- array with categories
    
    Returns:
        m -- line slope (y = m*x +b)
        b -- line intercept of the line
        std_error -- standard error of regression
        r_squared -- coefficient of determination
    
    Example:
        import numpy as np
        x = np.random.rand(10000)
        y = np.random.rand(10000)
        slope,intercept,std_error,r_squared = stochastic_theil_sen(x,y,200,20)
    """

    num_points = arr_x.shape[0]
    
    sampled_x = np.zeros(sample_size)
    sampled_y = np.zeros(sample_size)
    arr_slope = np.zeros(num_samples)
    arr_intercept = np.zeros(num_samples)
    std_error = 0.
    r_squared = 0.
        

    for i in xrange(num_samples):        
        population = np.random.randint(0, num_points, size = sample_size)
        for j in xrange(sample_size):
            sampled_x[j] = arr_x[population[j]]
            sampled_y[j] = arr_y[population[j]]

        med,b,s,r = theil_sen(sampled_x, sampled_y)
        std_error += s
        r_squared += r
        print i, med, b
        arr_slope[i] = med
        arr_intercept[i] = b
    
    slope = np.median(arr_slope)
    intercept = np.median(arr_intercept)
    std_error = std_error / num_samples
    r_squared = r_squared / num_samples
    
    return slope, intercept, std_error, r_squared


def eliminate_noise(slope, intercept, arr_x, arr_y, treshold=0.05):
    """
    Filter all data points that are too far from the line (treshold) defined by linear coefficients (slope, intercept).
    
    Arguments:
        slope -- line slope (y = slope * x + intercept)
        intercept -- line intercept
        arr_x -- array with x coordinates
        arr_y -- array with y coordinates
        treshold -- defines what is too far from the line (percentage from the line)
    
    Returns:
        arr_x -- array with x coordinates (filtered dataset)
        arr_y -- array with y coordinates (filtered dataset)
    
    Example:
        import numpy as np
        x = np.random.rand(10000)
        y = np.random.rand(10000)
        slope,intercept,std_error,r_squared = stochastic_theil_sen(x,y,200,20)
        x, y = eliminate_noise(slope, intercept, x, y, 0.1)
    """
    lst_x = []
    lst_y = []
    
    initial_size = arr_x.shape[0]
    
    def error(i):
        return np.absolute(arr_y[i] - (slope * arr_x[i] + intercept)) / arr_y[i]
 
    for i in xrange(initial_size):
        if error(i) < treshold:
            lst_x.append(arr_x[i])
            lst_y.append(arr_y[i])
    
    arr_x = np.asarray(lst_x)
    arr_y = np.asarray(lst_y)

    print '[Eliminate noise] Number of points: %d, shrinked to %.2f%% of %d' % (arr_x.shape[0], (arr_x.shape[0]*100.)/initial_size, initial_size)    

    return arr_x, arr_y

def cart_2_spher(x,y,z):
    """Convert the Cartesian vector [x, y, z] to spherical coordinates [r, theta, phi].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.
    """

    vector = [x,y,z]

    # The radial distance.
    r = norm(vector)

    # Unit vector.
    unit = vector / r

    # The polar angle.
    theta = m.acos(unit[2])

    # The azimuth.
    phi = m.atan2(unit[1], unit[0])

    # Return the spherical coordinate vector.
    return r, theta, phi


def spher_2_cart(r, theta, phi):
    """Convert the spherical coordinate vector [r, theta, phi] to the Cartesian vector [x, y, z].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.
    """

    # Trig alias.
    sin_theta = m.sin(theta)

    # The vector.
    x = r * m.cos(phi) * sin_theta
    y = r * m.sin(phi) * sin_theta
    z = r * m.cos(theta)

    return x, y, z


"""
def cart_2_spher(x, y, z):
    
    # transforms carthesian coordinates x, y, z
    # to spherical r, theta, phi    
            
    r = m.sqrt(x**2 + y**2 + z**2)
    theta = m.atan(m.sqrt(x**2 + y**2) / z)
    phi = m.atan(y/x)    

    return r, theta, phi
    
def spher_2_cart(r, theta, phi):
    
    # transforms spherical coordinates r, theta, phi    
    # to carthesian x, y, z
        
    x = r * m.sin(theta) * m.cos(phi)
    y = r * m.sin(theta) * m.sin(phi)
    z = r * m.cos(theta)
  
    return x, y, z
"""

def pixel_offset(row, column, num_stripes, pixel_width):
    """
    calculates pixel offset (dx,dy) from the center of the detector
    pixels are enumerated from 0 to n-1
    top left pixel is (0,0)
    """
    center = num_stripes / 2.    
    dx = (center - column - 0.5) * pixel_width
    dy = (center - row - 0.5) * pixel_width
    
    return dx, dy


def pixel_coordinates(row, column, num_stripes, pixel_width, detector_center, theta_d, phi_d):
    """
    calculate spherical coordinates for a pixel
    pixels are enumerated from 0 to n-1
    top left pixel is (0,0)
    front adcs are in row
    back adcs are in columns
    """
    dx, dy = pixel_offset(row, column, num_stripes, pixel_width)
    x = detector_center[0] + dx * m.cos(theta_d) 
    y = detector_center[1] + dy
    z = detector_center[2] - dx * m.sin(theta_d) * m.cos(phi_d)
    r, theta, phi = cart_2_spher(x, y, z)

    cartesian = [x, y, z]
    spherical = [r, theta, phi]
    return cartesian, spherical

def two_body_kinematics(theta, E_p, m_p, m_t, m_1, m_2, Q, excitation_energies=[0.]):  
    """
    two_body_kinematics
    
    arguments:
    theta - detector angle
    E_p - projectile energy
    m_p - projectile mass
    m_t - target mass
    m_1 - particle 1 mass
    m_2 - particle 2 mass
    Q - Q value, Q = 0 means elastic scattering, 
    excitation_energies = list of excitation energies, default is 0
    
    results:
    E_1 - particle 1 energy
    E_2 - particle 2 energy
    """     
     
    # helper variables
    alpha = ((m_2 - m_p) * E_p + m_2 * Q) / (m_1 + m_2)
    beta = m.sqrt(m_p * m_1 * E_p) / (m_1 + m_2)
    gamma = m.cos(theta)
    delta = m.sqrt(gamma**2 + alpha / beta**2)
        
    E_1 = alpha + 2 * beta**2 * gamma * (gamma + delta)
    E_2 = E_p + Q - E_1
    
    theta_2 = np.arccos((m_1 / (2.0 * np.sqrt(m_p * m_2 * E_p * E_2))) * (m_p * E_p / m_1 + m_2 *E_2 /m_1 - E_1))
    
    
    
    return E_1, E_2, theta_2


@mem.cache
def pixel_energies(detector, kinematics, target, stragg_projectile, stragg_target):
    '''
    izracunavamo koliku energiju ima cestica detektirana u nekom pixelu
    
    detector: num_stripes, detector_width, detector_distance, theta, phi
    theta je kut centra detektora, unosi se u stupnjevima
    phi se u ovakovom postavu uzima uvijek 0, s time da se onda theta za treci i cetvrti detektor uzima s predznakom minus
    broj stripova je i za horizintalne i vertikalne
    detector_width je sirina detektora dana u metrima
    detector_distance je udaljenost detektora od mete u metrima

    kinematics: E_p, m_p, m_t, m_1, m_2    
    parametri koji nam ulaze u kinematiku reakcije: 
    energija projektila (umanjena za gubitke do pola mete), masa projektila i mete, masa izlaznih cestica

    target: element, rmax, density
    parametri koji nam ulaze u stragg:
    element se odnosi na metu, ako se sastoji samo od jednog elementa
    rmax je udaljenost do pola mete dana u mg/cm2 - jer cestica nastala u reakciji mora
    proci jos pola mete/cos(theta_pixela) da dodje do odgovarajuceg pixela

    gustoca je defaultna iz stragga:
    za zlato gustoca je 19.3
    za bor gustoca je  2.350
    '''


    num_stripes = detector['num_stripes']
    detector_width = detector['detector_width']
    detector_distance = detector['detector_distance']
    theta = detector['theta']
    phi = detector['phi']
    
    E_p = kinematics['E_p']
    m_p = kinematics['m_p']
    m_t = kinematics['m_t']
    m_1 = kinematics['m_1']
    m_2 = kinematics['m_2']

    rmax = target['rmax']
    density = target['density']
    

    theta_r = m.radians(theta)
    phi_r = m.radians(phi)
    
    pixel_width = detector_width / num_stripes
    
    detector_center = spher_2_cart(detector_distance, theta_r, phi_r)

    final_energies = np.zeros([num_stripes,num_stripes])
    
    for row in xrange(num_stripes):
        for column in xrange(num_stripes):
    
            cartesian, spherical = pixel_coordinates(row, column, num_stripes, pixel_width, detector_center, theta_r, phi_r)
            pixel_theta = spherical[1]      
            
            E_1, E_2 = two_body_kinematics(pixel_theta, E_p, m_p, m_t, m_1, m_2)
            
            #print pixel_theta, E_p, m_p, m_t, m_1, m_2, E_1, E_2
            
            tmp_projectile = stragg_projectile + ', %f'%(E_1)
    
            # calling stragg
            dist = rmax/m.cos(pixel_theta)
            final_energies[row][column] = calculate_stragg_energy(tmp_projectile, stragg_target, dist, density)            
            
            
            print row*16+(column+1), row, column, E_1, final_energies[row][column], dist, tmp_projectile, stragg_target, density
            
            #print "|",
            
            """
            # mystragg num_elements e1 theta file symbol a2 rmax density
            proc = subprocess.Popen(
                                    ['./mystragg',
                                    '%f'%(E_1),
                                    '%s'%(element),
                                    '%d'%(m_t),
                                    '%f'%(rmax/m.cos(pixel_theta)),
                                    '%f'%(density)],
                                    stdout=subprocess.PIPE,shell=False)
            (out, err) = proc.communicate() # collect output
            matches = regex.search(out) # search for the final energy number
            final_energies[row][column] = float(matches.group(1)) # get the number and convert it to float

            #print matches.group(1)
            """

    print "\npixel_energies finished"
    
    return final_energies

def calculate_stragg_energy(projectile, target, rmax, density):
    """
    use stragg to calculate final energy
    
    Arguments:
    projectile = string - form 'symbol', A, Ein
    target = [stragg_string] - form: 'symbol', A
    target = [string1, string2, ..., string_n] - form: 'symbol', A, ratio
    rmax
    density
    
    Result:
    final_energy
    """

    child = pexpect.spawn('./stragg2')
    
    # Type in Projectile symbol 'quoted', A, Ein(MeV)
    child.expect('Type in Projectile',2)
    child.sendline(projectile)
    
    # How many elements in the target?
    child.expect('How many elements',2)
    num_elements = len(target)
    answer = "%d"%(num_elements)
    child.sendline(answer)
    
    # Enter 'symbol', and A of target
    if num_elements == 1:
        child.expect('Enter',2)
        child.sendline(target[0])
    else:
        for element in target:
            child.expect('Enter',2)
            child.sendline(element)
    
    # Enter target thickness (+ for mg/cm2, - for microns)
    child.expect('Enter target',2)
    child.sendline(str(rmax))
    
    # Enter target density, g/cm3 (default= 19.300)
    child.expect('Enter target density',2)
    child.sendline(str(density))
    
    # Full printout (Y/N)?
    child.expect('Full printout',2)
    child.sendline('N')
    
    """
     Slowing down of B (Z =    5 A =  10.00) in Au(Z =   79 A =  197.00)
    
     Incident energy =    49.00000 MeV
     Target density  =    19.30000 g/cm3
    
    
     Final energy=                     33.249 MeV
     Total range traversed=         19300.000 ug/cm2
                          =            10.000 um
     Energy loss=                      15.751 MeV
     Ionization energy deposited=      15.742 MeV
     L-S  straggling  width =    170.248 keV
     Corr. straggling width =   1027.789 keV
     total straggling width =   1041.794 keV
     ***********************************************
     new energy = (0=STOP)
    """
    child.expect('new energy',2)
    result = child.before
    child.close(force=True)
    
    regex = re.compile(r'.*?Final energy=\s+([\d.]+)')
    matches = regex.search(result) # search for the final energy number
    
    return float(matches.group(1)) # get the number and convert it to float


def calculate_cos_theta_12(theta_1, theta_2, phi_1, phi_2):
    """
    calculate relative angle between two detected particles
    """
    return np.cos(theta_1) * np.cos(theta_2) + np.sin(theta_1) * np.sin(theta_2) * np.cos(phi_1 - phi_2)


def get_excitation_energy(m_p,m_1,m_2,m_3,E_1,E_2,theta_1,theta_2,phi_1,phi_2,E_threshold):
    """
    calculate excitation energies
    kinemat equation 122
    """    
    tmp_0 = E_1/m_1 + E_2/m_2
    tmp_1 = np.sqrt(E_1 * E_2 / (m_1 * m_2))
    
    tmp_2 = calculate_cos_theta_12(theta_1, theta_2, phi_1, phi_2)
    tmp_3 = 1. * m_1 * m_2 / (m_1 + m_2)
    E_rel = tmp_3 * (tmp_0 - 2 * tmp_1 * tmp_2)
    
    return E_rel + E_threshold      
        

def get_excitation_energies(l_mass, l_energy, l_theta, l_phi, l_energy_thr):
    """
    calculate excitation energies
    
    Arguments:
    l_mass = [m1, m2, m3, mp]
    l_energy = [E1, E2, E3, Ep]
    l_theta = [theta1, theta2, theta3, thetap]
    l_phi [phi1, phi2]
    l_energy_thr = [Ethr12, Ethr13, Ethr23] # treshold energies
    
    Result:
    l_excitation = [Ex1, Ex2, Ex3]
    """
    def __tmp_mass(index):    
        tmp_mass_a = l_mass[index] / (l_mass[2] * (l_mass[index] + l_mass[2]))
        tmp_mass_b = (l_mass[index] + l_mass[2]) / l_mass[index]
        return [tmp_mass_a, tmp_mass_b]
    
    
    mu01 = 1. * l_mass[0] * l_mass[1] / (l_mass[0] + l_mass[1])
    cos_theta01 = np.cos(l_theta[0]) * np.cos(l_theta[1]) + np.sin(l_theta[0]) * np.sin(l_theta[1]) * np.cos(l_phi[0] - l_phi[1])
    
    E01 = mu01 * (l_energy[0] / l_mass[0] + l_energy[1] / l_mass[1] 
        - 2 * np.sqrt((l_energy[0] * l_energy[1]) / (l_mass[0] * l_mass[1]))
        * cos_theta01)

    tmp_mass_0a, tmp_mass_0b = __tmp_mass(0)
    tmp_mass_1a, tmp_mass_1b = __tmp_mass(1)
    
    tmp_sqrt0 = np.sqrt(l_mass[0] * l_energy[0])
    tmp_sqrt1 = np.sqrt(l_mass[1] * l_energy[1])
    tmp_sqrtp = np.sqrt(l_mass[3] * l_energy[3])

    E02 = tmp_mass_0a * ((tmp_mass_0b * tmp_sqrt0 * np.cos(l_theta[0])
            + tmp_sqrt1 * np.cos(l_theta[1]) - tmp_sqrtp)**2 
            + (-1 * tmp_mass_0b * tmp_sqrt0 * np.sin(l_theta[0])
            + tmp_sqrt1 * np.sin(l_theta[1]))**2)

    E12 = tmp_mass_1a * ((tmp_mass_1b * tmp_sqrt1 * np.cos(l_theta[1])
            + tmp_sqrt0 * np.cos(l_theta[0]) - tmp_sqrtp)**2 
            + (-1 * tmp_mass_1b * tmp_sqrt1 * np.sin(l_theta[1])
            + tmp_sqrt0 * np.sin(l_theta[0]))**2)
    
    result = [E01 + l_energy_thr[0],
            E02 + l_energy_thr[1], 
            E12 + l_energy_thr[2]]    
    
    return result

"""
Solid Angle of a Cartesian Surface Element
input:
  row - row index 0 ... num_stripes - 1
  column - column index 0 ... num_stripes - 1
  detector_width - width of the detector
  num_stripes - number of rows and columns
  R - detector distance
output:
  solid_angle - solid angle of the pixel defined by row and column
notes:
  http://web.utk.edu/~rpevey/NE406/lesson2.htm
"""
def calculate_solid_angle(row, column, detector_width, num_stripes, R):

    def __integral(W, L, z):
        temp = z * np.sqrt(W*W + z*z + L*L)
        temp = L * W / temp
        temp = np.arctan(temp)
        return temp 
        
    x1 = (row - num_stripes/2.) * detector_width/num_stripes 
    x2 = (row - num_stripes/2. + 1) * detector_width/num_stripes 
    y1 = (column - num_stripes/2.) * detector_width/num_stripes 
    y2 = (column - num_stripes/2. + 1) * detector_width/num_stripes 
    solid_angle = __integral(x2,y2,R) - __integral(x2,y1,R) - __integral(x1,y2,R) + __integral(x1,y1,R)
    
    return solid_angle


"""
Calculating 8Be energy and impulse from 2 alpha particles
angles are in degrees
input:
alpha_1 ... [energy, theta, phi]
alpha_2 ... [energy, theta, phi]
e_12
output:
8Be [energy, theta, phi]
"""
def calculate_8Be(alpha_1, alpha_2, E_DEC):
    
        
    E_8Be = alpha_1[0] + alpha_2[0] - E_DEC
    c_1 = np.sqrt(alpha_1[0]) * np.cos(np.deg2rad(alpha_1[1])) + np.sqrt(alpha_2[0]) * np.cos(np.deg2rad(alpha_2[1]))
    c_theta_8Be = c_1 / np.sqrt(2 * E_8Be)
    
    r_theta_8Be = np.arccos(c_theta_8Be)
    d_theta_8Be = np.rad2deg(r_theta_8Be)
    c_21 = np.sqrt(alpha_1[0]) * np.sin(np.deg2rad(alpha_1[1])) * np.sin(np.deg2rad(alpha_1[2]))
    c_22 = np.sqrt(alpha_2[0]) * np.sin(np.deg2rad(alpha_2[1])) * np.sin(np.deg2rad(alpha_2[2]))
    
    s_phi_8Be = (c_21 + c_22) / (np.sin(r_theta_8Be) * np.sqrt(2 * E_8Be))
    
    r_phi_8Be = np.arcsin(s_phi_8Be)
    d_phi_8Be = np.rad2deg(r_phi_8Be)
    
    return [E_8Be, d_theta_8Be, d_phi_8Be]




class TreeGenerator:    
    
    """ container for out_t arrays """
    def __init__(self, max_particles, file_new_tree):
        self.max_particles = max_particles        
        
        # zapisat cemo koincidencije na nacin da se sve cestice nadu u istom zapisu
        # poredane na nacin na koji su poredani detektori u LST_DET
        self.event = np.zeros(1, dtype='int')# broj eventa, varijabla tipa int32
        self.detector = np.zeros(max_particles, dtype='short')# ukljuceni detektori
        self.cnc = np.zeros(1, dtype='short') # oznaka koincidencija 1,2,3,4...
        self.wm = np.zeros(1, dtype='short') # multiplicitet eventa, varijabla tipa Short
        self.adc = np.zeros(max_particles*3, dtype='short') # oznaka piksela: front, back, delta_E
        self.ampl = np.zeros(max_particles*3, dtype='short') # amplituda
        self.nrg = np.zeros(max_particles*3, dtype='float') # energy levels
        
        self.ptype = np.zeros(max_particles, dtype='short') # particle type
        self.cnrg = np.zeros(max_particles, dtype='float') # corrected energy levels
        
        
        self.radius = np.zeros(max_particles, dtype='float') # angles for particles
        self.theta = np.zeros(max_particles, dtype='float') # angles for particles
        self.phi = np.zeros(max_particles, dtype='float') # angles for particles
        
        self.tree = self.__create_out_tree(file_new_tree)


    def __create_out_tree(self, file_new_tree):
        """ tree initialization """
        
        tree = TTree('T', '%s'%(file_new_tree)) # novo stablo
        
        # kreiramo strukturu stabla da bude identicna pocetnoj
        tree.Branch('event', self.event, 'event/l') # grana koja sadrzi multiplicitet
        tree.Branch('cnc', self.cnc, 'cnc/S') # oznaka koincidencija
        tree.Branch('wm', self.wm, 'wm/S') # grana koja sadrzi multiplicitet
        
        if self.max_particles == 1:        
            tree.Branch('detector', self.detector, 'detector/S') # oznaka detektora
            tree.Branch('radius', self.radius, 'radius[1]/D') # polje kutova
            tree.Branch('theta', self.theta, 'theta[1]/D') # polje kutova
            tree.Branch('phi', self.phi, 'phi[1]/D') # polje kutova
            tree.Branch('ptype', self.ptype, 'ptype[1]/S') # oznaka vrsta cestica
            tree.Branch('cnrg', self.cnrg, 'cnrg[1]/D') # polje korektiranih energija
        else:
            tree.Branch('detector', self.detector, 'detector[cnc]/S') # oznaka detektora
            tree.Branch('radius', self.radius, 'radius[cnc]/D') # polje kutova
            tree.Branch('theta', self.theta, 'theta[cnc]/D') # polje kutova
            tree.Branch('phi', self.phi, 'phi[cnc]/D') # polje kutova
            tree.Branch('ptype', self.ptype, 'ptype[cnc]/S') # oznaka vrsta cestica
            tree.Branch('cnrg', self.cnrg, 'cnrg[cnc]/D') # polje korektiranih energija
            
        tree.Branch('wadc', self.adc, 'wadc[wm]/S') # polje adc-a duljine wm
        tree.Branch('wampl', self.ampl, 'wampl[wm]/S') # polje amplituda duljine wm
        tree.Branch('wnrg', self.nrg, 'wnrg[wm]/D') # polje energija
        

    
        return tree
        
    def fill_event(self, oldtree):
        """ fill the tree """
        # kopiranje
        self.event[0] = oldtree.event
        self.cnc[0] = oldtree.cnc
        self.wm[0] = oldtree.wm
        
        if self.max_particles==1:
            self.detector[0] = oldtree.detector
            self.radius[0] = oldtree.radius
            self.theta[0] = oldtree.theta
	    self.phi[0] = oldtree.phi
	    try:
	    	self.ptype[0] = oldtree.ptype
	    except:
		pass
	    try: 
	    	self.cnrg[0] = oldtree.cnrg
	    except:
		pass
	    #self.ptype[0] = oldtree.ptype
            #self.cnrg[0] = oldtree.cnrg
        else:        
            for j in range(oldtree.cnc):
                self.detector[j] = oldtree.detector[j]                
                self.radius[j] = oldtree.radius[j]
                self.theta[j] = oldtree.theta[j]
                self.phi[j] = oldtree.phi[j]
		try:
                	self.ptype[j] = oldtree.ptype[j]
		except:
			pass
		try:
                	self.cnrg[j] = oldtree.cnrg[j]
		except:
			pass

        for j in range(oldtree.wm):
            self.adc[j] = oldtree.wadc[j]
            self.ampl[j] = oldtree.wampl[j]
            self.nrg[j] = oldtree.wnrg[j]

        #self.tree.Fill()
        
    def save_new_event(self):
        self.tree.Fill()
    
    def get_tree(self):
        """ fetch the tree """
        return self.tree
    
    
