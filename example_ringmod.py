from technologies import silicon_photonics
from ringmodulator import AddDropRingModulator
from picazzo3.traces.wire_wg.trace import WireWaveguideTemplate

wg_tmpl = WireWaveguideTemplate()

ad = AddDropRingModulator(trace_template=wg_tmpl)
ad_lay = ad.Layout()

ad_lay.visualize()
