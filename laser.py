import ipkiss3.all as i3
from numpy import real, exp, sqrt, log
from scipy.constants import speed_of_light

class Laser(i3.PCell):

    class Netlist(i3.NetlistView):
        # To build a netlist, check out the netlist guide: http://docs.lucedaphotonics.com/3.1/guides/netlist/index.html
        def _generate_terms(self, terms):
            terms += i3.ElectricalTerm(name='in')
            terms += i3.OpticalTerm(name='out')
            return terms

    class CapheModel(i3.CapheTODEModelView):
        """Single-mode laser diode model.

        This node describes a laser diode. It has two ports, the first one is the
        electrical current input (in A), the second one is the optical output.

        The lasing wavelength is assumed to be the same as the system wavelength.
        The laser itself can induce a frequency shift, this is translated in 
        a rapid varying phase change in the optical output.

        Note that typically the lasing threshold current is very temperature
        dependent. It should be calculated automatically in the attached model,
        i.e., the model is responsible of setting I_th based on the current laser 
        parameters.

        The default values were taken for a typical semiconductor laser diode,
        with a lasing threshold of about 17.8 mA.

        References
        ----------
        - Agrawal's book "Fiber-Optic Communication Systems" (3rd edition)
        - PhD of Geert De Pestel "Studie van halfgeleider laserdiodes onder analoge
          microgolfuitsturing"
        - Geert Morthier's invaluable input to help finding the correct parameters 
          and performing a few small fixes to the model.
        """

        beta = i3.FloatProperty(doc="Spontaneous emission factor (-)", default=5e-5)
        tau_e = i3.FloatProperty(doc="Electron lifetime (ns)", default=1.0)
        gamma = i3.FloatProperty(doc="Confinement factor (-)", default=0.35)
        v_g = i3.FloatProperty(doc="Group speed (m/s)", default=3e8/4.0)
        sigma_g = i3.FloatProperty(doc="Differential gain (cm^2)", default=2.5*1e-16)
        N_T = i3.FloatProperty(doc="Carrier density at transparency (cm^-3)", default=1.0e18)
        V = i3.FloatProperty(doc="Active volume (um^3)", default=0.1*2*300)
        epsilon_nl = i3.FloatProperty(doc="Gain compression factor (-)", default=0.0)
        n_sp = i3.FloatProperty(doc="Spontaneous emission factor (-)", default=2.0)
        
        # Laser losses, splitted into internal losses and mirror losses
        alpha_int = i3.FloatProperty(doc="Internal losses (including free-carrier absorption, scattering and other possible mechanisms) (cm^-1)", 
                                             default=25.0)
        alpha_mirror = i3.FloatProperty(doc="Mirror losses 1/(2L)*ln(1/R1*1/R2) (cm^-1)",
                                                default=1.0 / (2*300e-4) * log(1 / 1.0 * 1 / 0.05))
        
        # Linewidth enhancement factor (named after Charles Henry)
        alpha_h = i3.FloatProperty(doc="Linewidth enhancement factor (-)", default=4)
                                        
        # Calculated properties
        tau_p = i3.FloatProperty(doc="[calculated] Photon lifetime, 1/(v_g*alpha_cav), where alpha_cav = alpha_int + alpha_mirror (ps)")
        N_0 = i3.FloatProperty(doc="[calculated] Number of carriers in the active volume for transparency, equal to N_T*V")
        G_N = i3.FloatProperty(dohc="[calculated] Rate of stimulated emission per electron above transparency")
        I_th = i3.FloatProperty(doc="[calculated] Laser threshold, based on N_0, G_N, tau_e, tau_p, v_g and sigma_g")

        def _default_tau_p(self):
            # Calculate the photon lifetime based on the cavity losses
            # and mirror losses.        
            # For the definition of the units, check out the documentation of the properties. 
            # For example, alpha_cav is loss in cm-1. The photon lifetime is expressed in ps. Hence
            # the multiplication factors 1e2 and 1e12.
            alpha_cav = self.alpha_int + self.alpha_mirror
            return 1.0/(self.v_g * (alpha_cav*1e2))*1e12 

        def _default_N_0(self):
            return (self.N_T * 1e6) * (self.V * 1e-18)

        def _default_G_N(self):
            return self.gamma * self.v_g * (self.sigma_g * 1e-4) / (self.V * 1e-18)

        def _default_I_th(self):
            q = 1.602176565*1e-19
            N_th = self.N_0 + 1 / (self.G_N * (self.tau_p * 1e-12))
            return q * N_th / (self.tau_e * 1e-9)


        def _calculate_S(self, environment, term1, term2, mode1, mode2):
            return 0
        
        @i3.compile_function()
        def _calculate_signal_ext(self, environment, t, s_in, s_ext, y):
            h = 6.62606957*1e-34
            """        
            See Agrawal's book "Fiber-Optic Communication Systems" (3rd edition),
            chapter "Laser characteristics" -->  "CW Characteristics"
            for a derivation and side notes.
            i.e. for FP lasers with coated facets, or for DFB lasers, these
            equations have to be modified properly.
            
            (self.v_g * self.alpha_mirror) is the rate with which the
            photons escape from the two facets. We assume one of the
            facets has 100% reflection. If both facets would have a non 100%
            reflectivity, we'd need to split up in a left and right port.

            Due to carrier-induced changes in the mode index, the output phase
            changes. This is captured in the phase variable y[2]

            Note: wavelength is expressed in um. So we multiply with 1e-6 to get SI units.
            """
            
            S = real(y[1][0])
            phase = real(y[2][0])
            c = 299792458.0
            P_e = (self.v_g * self.alpha_mirror * 100) * h * 1.0 / (environment.wavelength * 1e-6) * S * c
            
            s_ext[0] = 0
            s_ext[1] = sqrt(P_e)*exp(1j*phase)

        
        @i3.compile_function()
        def _calculate_dydt(self, environment, t, y, dydt, s_in):
            """
            N = y[0], number of electrons
            S = y[1], number of photons

            s_in[0](t): Electrical current input in the laser diode
            s_in[1](t): Optical input (not used right now, the laser only emits
                        light but has no feedback into the laser)
            """
            q = 1.602176565*1e-19

            N = real(y[0][0])
            P = real(y[1][0])
            I = s_in[0][0]

            G_ = self.G_N * (N - self.N_0)
            G = G_ * (1 - self.epsilon_nl * P)

            """The spontaneous emission rate is tricky to calculate

            A good approximation is by using 
            R_sp = self.beta * N / (self.tau_e * 1e-9),
            but here beta is quite empyrically chosen and not really close
            to any physical meaning.

            Otherwise, it's more accurate to use 
            R_sp = self.n_sp * G,

            but then n_sp is a function of time, and it can even be a negative value
            (i.e., to counter the negative G when below lasing threshold).
            n_sp is then given by 1 / (1 - exp( (h*f - e*V) / (k*T) ))
            And e*V has to deal with the fermi levels,
            E_f_c - E_f_v = e*V (conduction band, valence band).

            When lasing, G is approximately equal to 1/tau_p. In this case,
            R_sp = 1/tau_p * n_sp

            This approximation could be used to ensure R_sp is always positive, 
            but it will not be accurate below lasing threshold.
            """
            #R_sp = self.n_sp * G # Use when you use noise
            R_sp = self.beta * N / (self.tau_e * 1e-9)

            # dN/dt
            dydt[0] = I / q - N / (self.tau_e * 1e-9) - G * P

            # dS/dt
            dydt[1] = G * P + R_sp - P / (self.tau_p * 1e-12)

            # d_phase/dt
            dydt[2] = 0.5 * self.alpha_h * (G_ - (1.0 / (self.tau_p * 1e-12)))


        def _calculate_variables(self, vt):
            vt.add_variable(id=2000, name="N, number of free electrons", nr_vars=1, scaling=1)
            vt.add_variable(id=2001, name="S, number of photons in the cavity mode", nr_vars=1, scaling=1)
            vt.add_variable(id=2002, name="phi, phase change due to carrier-induced changes in the mode index", nr_vars=1, scaling=1)
