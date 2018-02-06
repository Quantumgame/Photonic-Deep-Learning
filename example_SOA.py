# -*- coding: utf-8 -*-

#%%
from technologies import silicon_photonics

import ipkiss3.all as i3
from testbench import MySource, MyDetector, MyElectricalSource
from soa import SOA

import numpy as np
from pylab import plt

i3.RoundedWaveguide.set_default_view("CapheModelFromLayout")
i3.Waveguide.set_default_view("CapheModelFromLayout")

t0 = 0
dt = 2e-12
t1 = 40e-9
src_bitrate = 0.5e9

src = MySource()
#src.CapheModelOnOff(bitrate=src_bitrate)
src_cm = src.CapheModelPRBS(amplitude=0.03535, bitrate=src_bitrate, noise_amplitude=0.0)
src.set_default_view(src_cm)

el_src = MyElectricalSource()
el_src_cm = el_src.CapheModelOn(t_start=0, amplitude=0.18666)
el_src.set_default_view(el_src_cm)

det = MyDetector()

det.CapheModel(R=0.2, BW=10e9)

soa = SOA()
soa.CapheModel(tau_carrier=300e-12)

from picazzo3.routing.place_route import PlaceAndConnect

testbench = PlaceAndConnect(
    child_cells={
        'src': src,
        'el_src': el_src,
        'soa': soa,
        'det': det,
    },
    links=[('src:out', 'soa:in'), ('soa:out', 'det:in'), ('el_src:out', 'soa:electrical_in')]
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
plt.subplot(211)
plt.plot(times * 1e9, np.abs(sources[:, 0]), 'g', linewidth=1, label='Input')
plt.plot(times * 1e9, np.abs(outputs[:, 0]), 'r', linewidth=2, label='Output')
plt.legend()
plt.subplot(212)
plt.plot(times * 1e9, np.real(states[:, 0]), 'b', label='Integrated gain SOA')
plt.legend()
plt.xlabel("Time ($\mu$ s)")
plt.show()
