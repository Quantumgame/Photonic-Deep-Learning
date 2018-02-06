"""Implementation of a deep neural network.
"""
 
#from picazzo3.routing.place_route import PlaceAndAutoRoute
from auto_place_and_connect import AutoPlaceAndConnect
from ringmodulator import AddDropRingModulator
import ipkiss3.all as i3

class DeepNeuralNetwork(AutoPlaceAndConnect):
    input_array = i3.ChildCellProperty(doc="The input array")
    rings = i3.ChildCellListProperty(doc="List of ring resonators")
    carrier_wgs = i3.ChildCellListProperty(doc="Bus waveguide pieces between the rings")
    source_wgs = i3.ChildCellListProperty(doc="Source waveguide pieces between the rings")

    nr_rows = i3.PositiveIntProperty(locked=True, doc="Number of rows in the input pixel array")
    nr_cols = i3.PositiveIntProperty(locked=True, doc="Number of rows in the input pixel array")
    nr_rings = i3.PositiveIntProperty(locked=True, doc="Number of ring resonators, equal to nr_rows * nr_cols")

    trace_template = i3.TraceTemplateProperty()

    def _default_nr_rings(self):
        return self.input_array.nr_rows * self.input_array.nr_cols

    def _default_nr_rows(self):
        return self.input_array.nr_rows

    def _default_nr_cols(self):
        return self.input_array.nr_cols

    def _default_rings(self):
        return [AddDropRingModulator(trace_template=self.trace_template) for x in xrange(self.nr_rings)]

    def _default_source_wgs(self):
        return [i3.Waveguide(trace_template=self.trace_template) for x in xrange(self.nr_rings)]

    def _default_carrier_wgs(self):
        return [i3.Waveguide(trace_template=self.trace_template) for x in xrange(self.nr_rings)]

    def _default_child_cells(self):
        # Define your child cells here. For example:
        child_cells = {'input_array': self.input_array}
        for x in xrange(1, self.nr_rows + 1):
            for y in xrange(1, self.nr_cols + 1):
                nr = (x - 1) * self.nr_cols + y - 1
                child_cells.update({'ring_{}_{}'.format(x, y): self.rings[nr]})

                if x < self.nr_rows and y < self.nr_cols:
                    a = 1 if y == self.nr_cols else 0
                    b = 1 if y < self.nr_cols else 0
                    wg_index = "{}{}_{}{}".format(x, y, x + a, y + b)

                    child_cells.update({'wg_carrier_{}'.format(wg_index): self.carrier_wgs[nr]})
                    child_cells.update({'wg_source_{}'.format(wg_index): self.source_wgs[nr]})
        return child_cells

    def _default_links(self):
        # Define the links between the components here as a list of tuples.
        # For example:
        links = []
        for x in xrange(1, self.nr_rows + 1):
            for y in xrange(1, self.nr_cols + 1):
                nr = (x - 1) * self.nr_cols + y - 1
                if x < self.nr_rows and y < self.nr_cols:
                    a = 1 if y == self.nr_cols else 0
                    b = 1 if y < self.nr_cols else 0

                    ring1 = "ring_{}_{}".format(x, y)
                    ring2 = "ring_{}_{}".format(x + a, y + b)
                    wg_index = "{}{}_{}{}".format(x, y, x + a, y + b)

                    links.append(('{}:out1'.format(ring1), 'wg_source_{}:in'.format(wg_index)))
                    links.append(('{}:in1'.format(ring2), 'wg_source_{}:out'.format(wg_index)))

                    links.append(('{}:in2'.format(ring1), 'wg_carrier_{}:in'.format(wg_index)))
                    links.append(('{}:out2'.format(ring2), 'wg_carrier_{}:out'.format(wg_index)))

        return links

    class Layout(AutoPlaceAndConnect.Layout):

        def _default_carrier_wgs(self):
            wgs = []
            for i, t in enumerate(self.cell.carrier_wgs):
                t = t.get_default_view(i3.LayoutView)
                t.set(shape=[(0, 0), (0, 10)])
                wgs.append(t)
            return wgs

        def _default_source_wgs(self):
            wgs = []
            for i, t in enumerate(self.cell.source_wgs):
                t = t.get_default_view(i3.LayoutView)
                t.set(shape=[(0, 0), (0, 10)])
                wgs.append(t)
            return wgs

        def _default_child_transformations(self):
            # Define the child transformations here
            return {'input_array': (0, 0)}

