from numpy.testing import assert_array_equal


class TestActiveLinkGradients(object):
    def setup_grids(self):
        self.grids = {
            'unit_spacing': RasterModelGrid(4, 5),
            'non_unit_spacing': RasterModelGrid(4, 5, 2.),
        }

    def test_unit_spacing(self):
        (rmg, values_at_nodes) = RasterModelGrid(4, 5), np.arange(20.)
        grads = rmg.calculate_gradients_at_active_links(values)

        assert_array_equal(grads, np.array([5, 5, 5, 5, 5, 5, 5, 5, 5,
                                            1, 1, 1, 1, 1, 1, 1, 1]))

        diffs = rmg.calculate_diff_at_active_links(values)
        assert_array_equal(grads, diffs)
        
    def test_non_unit_spacing(self):
        _SPACING = 2.
        (rmg, values_at_nodes) = RasterModelGrid(4, 5, _SPACING), np.arange(20.)
        grads = rmg.calculate_gradients_at_active_links(values)
        assert_array_equal(grads, (1. / _SPACING) *
                              np.array([5, 5, 5, 5, 5, 5, 5, 5, 5,
                                        1, 1, 1, 1, 1, 1, 1, 1]))
        diffs = rmg.calculate_diff_at_active_links(values)
        assert_array_equal(grads, (1. / _SPACING) * diffs)

    def test_out_array(self):
        (rmg, values_at_nodes) = RasterModelGrid(4, 5), np.arange(20.)
        grads = np.empty(17)
        rtn_grads = rmg.calculate_gradients_at_active_links(values,
                                                            out=grads)
        assert_array_equal(grads, np.array([5, 5, 5, 5, 5, 5, 5, 5, 5,
                                               1, 1, 1, 1, 1, 1, 1, 1]))
        self.assertIs(rtn_grads, grads)
        
    def test_diff_out_array(self):
        rmg = RasterModelGrid(4, 5)
        values = np.arange(20)
        diff = np.empty(17)
        rtn_diff = rmg.calculate_diff_at_active_links(values, out=diff)
        assert_array_equal(
            diff,
            np.array([5, 5, 5, 5, 5, 5, 5, 5, 5,
                      1, 1, 1, 1, 1, 1, 1, 1]))
        self.assertIs(rtn_diff, diff)


class TestLinkGradients(unittest.TestCase, NumpyArrayTestingMixIn):
    def test_unit_spacing(self):
        rmg = RasterModelGrid(4, 5)
        values = np.arange(20, dtype=float)
        grads = rmg.calculate_gradients_at_links(values)
        assert_array_equal(
            grads,
            np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     dtype=float))
        diffs = rmg.calculate_diff_at_links(values)
        assert_array_equal(grads, diffs)
        
    def test_non_unit_spacing(self):
        rmg = RasterModelGrid(4, 5, 2.)
        values = np.arange(20, dtype=float)
        grads = rmg.calculate_gradients_at_links(values)
        assert_array_equal(
            grads,
            .5 * np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                          dtype=float))
        diffs = rmg.calculate_diff_at_links(values)
        assert_array_equal(grads, .5 * diffs)

    def test_out_array(self):
        rmg = RasterModelGrid(4, 5)
        values = np.arange(20, dtype=float)
        grads = np.empty(31.)
        rtn_grads = rmg.calculate_gradients_at_links(values, out=grads)
        assert_array_equal(
            grads,
            np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     dtype=float))
        self.assertIs(rtn_grads, grads)

    def test_diff_out_array(self):
        rmg = RasterModelGrid(4, 5)
        values = np.arange(20, dtype=float)
        diff = np.empty(31.)
        rtn_diff = rmg.calculate_diff_at_links(values, out=diff)
        assert_array_equal(
            diff,
            np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     dtype=float))
        self.assertIs(rtn_diff, diff)
        
        
class TestMaxGradients(unittest.TestCase, NumpyArrayTestingMixIn):
    def test_scalar_arg(self):
        rmg = RasterModelGrid(4, 5)
        values = np.arange(20, dtype=float)
        grad = rfuncs.calculate_max_gradient_across_cell_faces(rmg, values, 0)
        self.assertEqual(grad, 5.)

    def test_iterable(self):
        rmg = RasterModelGrid(4, 5)
        values = np.arange(20, dtype=float)
        grad = rfuncs.calculate_max_gradient_across_cell_faces(rmg, values, [0, 4])
        assert_array_equal(grad, [5., 5.])

    def test_scalar_arg_with_links(self):
        rmg = RasterModelGrid(4, 5)
        values = np.array([0, 1,  3, 6, 10,
                           0, 1,  3, 6, 10,
                           0, 1,  3, 5, 10,
                           0, 1, -3, 6, 10,])
        (grad, face) = rfuncs.calculate_max_gradient_across_cell_faces(
            rmg, values, (0, 4), return_face=True)
        assert_array_equal(grad, [1, 6])
        assert_array_equal(face, [1, 2])

        link_ids = rfuncs.active_link_id_of_cell_neighbor(rmg, [0, 4], face)
        assert_array_equal(link_ids, [9, 7])

        node_ids = rfuncs.node_id_of_cell_neighbor(rmg, [0, 4], face)
        assert_array_equal(node_ids, [5, 17])

class TestGradientsAcrossFaces(unittest.TestCase, NumpyArrayTestingMixIn):
    def test_no_arg(self):
        rmg = RasterModelGrid(4, 5)
        values = np.array([0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,])
        grads = rfuncs.calculate_gradients_across_cell_faces(rmg, values)
        assert_array_equal(
            grads,
            np.array([[-2, 0, 1, 0], [-3, 0, 2, 0], [-4, 0, 3, 0],
                      [-2, 0, 1, 0], [-3, 0, 2, 0], [-4, 0, 3, 0]]))

    def test_scalar_arg(self):
        rmg = RasterModelGrid(4, 5)
        values = np.array([0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,])
        grads = rfuncs.calculate_gradients_across_cell_faces(rmg, values, 0)
        assert_array_equal(grads, np.array([[-2, 0, 1, 0]]))

    def test_iterable_arg(self):
        rmg = RasterModelGrid(4, 5)
        values = np.array([0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,
                           0, 1, 3, 6, 10,])
        grads = rfuncs.calculate_gradients_across_cell_faces(
            rmg, values, np.array([0, 4]))
        assert_array_equal(
            grads, np.array([[-2, 0, 1, 0], [-3, 0, 2, 0]]))


class TestFluxDivergence(unittest.TestCase, NumpyArrayTestingMixIn):
    def test_flux_from_south_to_north(self):
        rmg = RasterModelGrid(4, 5)
        active_link_flux = np.array([0., 0., 0., 1., 1., 1., 3., 3., 3.,
                                     0., 0., 0., 0., 0., 0., 0., 0.])
        divs = rmg.calculate_flux_divergence_at_nodes(active_link_flux)

        assert_array_equal(
            divs,
            np.array([0.,  0.,  0.,  0., 0.,
                      0.,  1.,  1.,  1., 0.,
                      0.,  2.,  2.,  2., 0.,
                      0., -3., -3., -3., 0.]))
        
    def test_flux_from_east_to_west(self):
        rmg = RasterModelGrid(4, 5)
        active_link_flux = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                     0., 1., 3., 6., 0., 1., 3., 6.])
        divs = rmg.calculate_flux_divergence_at_nodes(active_link_flux)

        assert_array_equal(
            divs,
            np.array([0., 0., 0., 0.,  0.,
                      0., 1., 2., 3., -6.,
                      0., 1., 2., 3., -6.,
                      0., 0., 0., 0.,  0.]))
        
    def test_out_array(self):
        rmg = RasterModelGrid(4, 5)
        active_link_flux = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                     0., 1., 3., 6., 0., 1., 3., 6.])

        divs = np.empty(20)
        rtn_divs = rmg.calculate_flux_divergence_at_nodes(active_link_flux,
                                                          out=divs)

        assert_array_equal(
            divs,
            np.array([0., 0., 0., 0.,  0.,
                      0., 1., 2., 3., -6.,
                      0., 1., 2., 3., -6.,
                      0., 0., 0., 0.,  0.]))
        self.assertIs(rtn_divs, divs)


