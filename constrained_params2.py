from __future__ import division
from math import log

def constraint_vesc(logmchi, experiment, vesc = 756.0):
    #print logmchi
    mchi = 10**logmchi
    bits = experiment.split("/")

    if 'XENON' in bits:
        # print 'Xenon'
        N_mass = (128-54)*0.939565560 + 54 * 0.938272046

    elif 'GERMANIUM' in bits:
        # print 'GERMANIUM'
	N_mass = (70-32)*0.939565560 + 32 * 0.938272046

    elif 'ARGON' in bits:
        # print 'Argon'
	N_mass = (40-18)*0.939565560 + 18 * 0.938272046

    else:
        N_mass =0
        print ('The constraint is not working')

    #print vesc
    c = 2.99792458e5
    x = mchi/N_mass
    reduced_mass = (N_mass* mchi)/(N_mass + mchi)
    lhs = (2.*(reduced_mass**2.)/N_mass)*((vesc/c))**2.*1.0e6

    return lhs
