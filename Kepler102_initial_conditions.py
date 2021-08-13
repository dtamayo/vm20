# run as python Kepler102_initial_conditions.py -s 100 -n 50
# where 100 is the index at which to start, i.e. first sim = run0100.bin
# 50 is the number of simulations to run, so last sim = run0149.bin

import numpy as np
import rebound
import astropy.units as u
import argparse
import pandas as pd

def _parse_args():
    parser = argparse.ArgumentParser('make_initial_conditions', 
                                     description='make_initial_conditions')

    parser.add_argument('-s', 
                        '--startindex', 
                        required=True, 
                        help='index at which to start')
    
    parser.add_argument('-n', 
                        '--nsims', 
                        required=True, 
                        help='number of simulations to generate initial conditions for')

    return parser.parse_args()

def Mmax(R, rhomax=15): # max mass in Mearth given R in Rearth and rhomax in g/cc
    V = 4/3*np.pi*(R*u.Rearth.to(u.cm))**3 # cm^3
    Mmax = 15*V*u.g.to(u.Mearth)
    return Mmax

def make_sim(idx, mb, mc, md, me, mf):
    mearth = 3e-6
    np.random.seed(idx)
    e = np.random.random()*0.08 + 0.01 # U(0.01, 0.09)
    sim = rebound.Simulation()
    sim.add(m=0.8)
    sim.add(m=mb[idx]*mearth, P=5.287, e=0, M=np.random.random()*2*np.pi)
    sim.add(m=mc[idx]*mearth, P=7.07, e=e, pomega=np.random.random()*2*np.pi, M=np.random.random()*2*np.pi)
    sim.add(m=md[idx]*mearth, P=10.31, e=e, pomega=np.random.random()*2*np.pi, M=np.random.random()*2*np.pi)
    sim.add(m=me[idx]*mearth, P=16.15, e=e, pomega=np.random.random()*2*np.pi, M=np.random.random()*2*np.pi)
    sim.add(m=mf[idx]*mearth, P=27.45, e=e, pomega=np.random.random()*2*np.pi, M=np.random.random()*2*np.pi)
    sim.move_to_com()
    sim.integrator = "whfast"
    sim.collision = "line"
    sim.dt = sim.particles[1].P/20
    sim.ri_whfast.safe_mode = 0
    sim.ri_whfast.keep_unsynchronized = 1
    for p in sim.particles[1:]:
        p.r = p.a * (p.m/sim.particles[0].m)**(1/3)
    sim.save('data/unfinished/run{0:04d}.bin'.format(idx))

def main():
    args = _parse_args()
    nsims = int(args.nsims)
    start = int(args.startindex)
    Rb = 0.5 # Rearth
    Rc = 0.62
    Rd = 1.31
    Re = 2.4
    Rf = 0.93
    Mmaxb = Mmax(Rb)
    Mmaxc = Mmax(Rc)
    Mmaxd = Mmax(Rd)
    Mmaxe = Mmax(Re)
    Mmaxf = Mmax(Rf)

    me = np.random.normal(loc=7.5, scale = 1.5, size=10000)

    mb = pd.read_csv('mb.txt', sep=' ')['x'].values
    mc = pd.read_csv('mc.txt', sep=' ')['x'].values
    md = pd.read_csv('md.txt', sep=' ')['x'].values
    mf = pd.read_csv('mf.txt', sep=' ')['x'].values

    mb = mb[mb < Mmaxb]
    mc = mc[mc < Mmaxc]
    md = md[md < Mmaxd]
    me = me[me < Mmaxe]
    mf = mf[mf < Mmaxf]

    for i in range(start, start + nsims):
        make_sim(i, mb, mc, md, me, mf)

if __name__ == '__main__':
    main()
