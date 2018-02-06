from technologies import silicon_photonics

from picazzo3.traces.wire_wg.trace import WireWaveguideTemplate
from DNN import DeepNeuralNetwork
from testbench import InputArray

wg_tmpl = WireWaveguideTemplate()

pixels = InputArray(nr_rows=3, nr_cols=3)

nn = DeepNeuralNetwork(trace_template=wg_tmpl,
                       input_array=pixels)

print nn.child_cells

nn_lay = nn.Layout()
print nn_lay.child_transformations
nn_lay.visualize()
