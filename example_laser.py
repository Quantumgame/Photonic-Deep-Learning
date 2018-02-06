# -*- coding: utf-8 -*-

#%%
from technologies import silicon_photonics

import ipkiss3.all as i3
from testbench import MyElectricalSource, MyDetector
from laser import Laser

import numpy as np
from pylab import plt

np.random.seed(101)

i3.RoundedWaveguide.set_default_view("CapheModelFromLayout")
i3.Waveguide.set_default_view("CapheModelFromLayout")

t0 = 0
dt = 1e-12
t1 = 20e-9
src_bitrate = 0.3e9

src = MyElectricalSource()
#src.CapheModelOnOff(bitrate=src_bitrate)
src_cm = src.CapheModelPRBS(amplitude=0.03535, bitrate=src_bitrate, noise_amplitude=0.0)
src.set_default_view(src_cm)

det = MyDetector()

det.CapheModel(R=0.2, BW=10e9)

laser = Laser()
laser_cm = laser.CapheModel(beta=5e-5,                       # Spontaneous emission factor (-)
                            tau_e=1.0,                       # Electron lifetime (ns)
                            gamma=0.35,                      # Confinement factor (-)
                            v_g=3e8/4.0,                     # Group speed (m/s)
                            sigma_g=2.5*1e-16,               # Differential gain (cm^2)
                            N_T=1.0e18,                      # Carrier density at transparency (cm^-3)
                            V=0.1 * 2 * 300,                 # Active volume (um^3)
                            epsilon_nl=0.0,                  # Gain compression factor (-) 
                            n_sp=2.0,                        # Spontaneous emission factor (not used for now)
                            alpha_int=25.0,                  # Internal losses (including free-carrier absorption, scattering and other possible mechanisms) (cm^-1)
                            alpha_mirror=1.0 / (2 * 300e-4)*np.log(1 / 1.0 * 1 / 0.05), # Mirror losses 1/(2L)*ln(1/R1*1/R2) (cm^-1)
                            alpha_h=4)                       # Linewidth enhancement factor (-)

print("Laser threshold: {} mA".format(laser_cm.I_th * 1e3)) # For the parameters given above, this should be 17.8450*1e-3

from picazzo3.routing.place_route import PlaceAndConnect

testbench = PlaceAndConnect(
    child_cells={
        'electrical_src': src,
        'laser': laser,
        'det': det,
    },
    links=[('electrical_src:out', 'laser:in'), ('laser:out', 'det:in')]
)
tb_cm = testbench.CapheModel()

# %%
# Simulation
#
from ipkiss3.simulation.engines.caphe_circuit_sim.caphenodegenerator import create_caphe_node
import caphe


env = caphe.base.EnvironmentObject(name='testenv', wavelength=1.55)
caphe_node = create_caphe_node(tb_cm, n_modes=1)


solver = caphe.base.CapheNodeSolver(caphe_node, env, check_scatter_matrices=False, check_linked=True, ignore_externals_not_linked=True, auto_flatten=True)

# Run simulation
print("Start simulation")
solver.set_integration_method(caphe.solvers.runge_kutta4)
solver.set_internal_dt(dt)
solver.solve(t0=t0, t1=t1, dt=dt, environment=env)
print("Done simulating")

# Retrieve results
times, states, outputs, sources = solver.get_states_and_output_and_sources()

#from output_mapping import get_outputs_map
#print get_outputs_map(caphe_node)
# %%
# Plotting
#
plt.subplot(311)
plt.plot(times * 1e9, np.abs(sources[:, 0]), 'g', linewidth=1, label='Electrical input')
plt.plot(times * 1e9, np.abs(outputs[:, 0]), 'r', linewidth=2, label='Output')
plt.legend()
plt.subplot(312)
plt.plot(times * 1e9, np.real(states[:, 0]), 'k', label='N (number of free electrons)')
plt.legend()
plt.subplot(313)
plt.plot(times * 1e9, np.real(states[:, 1]), 'k', label='S (number of photonics in the cavity mode)')
plt.legend()
plt.xlabel("Time (ns)")
plt.show()

#vt.add_variable(id=2000, name="N, number of free electrons", nr_vars=1, scaling=1)
#vt.add_variable(id=2001, name="S, number of photons in the cavity mode", nr_vars=1, scaling=1)
#vt.add_variable(id=2002, name="phi, phase change due to carrier-induced changes in the mode index", nr_vars=1, scaling=1)
