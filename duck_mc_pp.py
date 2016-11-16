# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 13:19:49 2016

@author: sthagon
"""

import argparse
import pp, sys, importlib



"""
hendlanje ulaznih argumenata
input:
    experiment module - python skripta koja vrti jednu iteraciju eksperimenta
    output_file - gdje se zapisuju rezultati
    detector_setup - oznaka postava detektora: postav_1, postav_2...
    number_of_workers - broj procesa koji ce raditi u paraleli
    e.g. python duck_mc_pp.py experiment_10B10B mc-results-10B10B-postav1.csv postav_1 2
"""

parser = argparse.ArgumentParser(description='Run Monte Carlo Experiments in Parallel')

parser.add_argument('experiment_module', type=str)

parser.add_argument('output_file', type=str)

parser.add_argument('detector_setup', type=str)

parser.add_argument('number_of_workers', default=1, type=int)

args = parser.parse_args()

print args


"""
import specificiranog experiment modula
"""
expmod = importlib.import_module(args.experiment_module)


"""
pokretanje i konfiguracija job servera
"""
job_server = pp.Server(secret='ZekoPeko&%!', ncpus=args.number_of_workers)
print "Starting pp with", job_server.get_ncpus(), "workers"
jobs = []


"""
MIJENJAMO ENERGIJU E_12 OD 0 DO E_MAX_REL
E_12 - energija iznad praga koji je potreban za raspad X -> C + D
"""
E_12 = expmod.E_12_START
while E_12 <= expmod.E_12_STOP:

    job = job_server.submit(
                        expmod.duck_mc_experiment, 
                        (E_12, expmod.BROJ_PONAVLJANJA, args.detector_setup), 
                        (), 
                        (expmod.DEPENDENCY_MODULE,))
    jobs.append(job)
    
    """
    povecavamo E_12 i prelazimo na sljedecu iteraciju izracuna
    """    
    E_12 += expmod.E_12_STEP


bufsize = 1
target = open(args.output_file, "w", bufsize)

# dohvacamo rezultate analiza
for job in jobs:
    result = job()
    print result
    target.write(result)
    target.write("\n")

# ispisujemo neke statistike
job_server.print_stats()


target.close()
