from scipy import weave


def erode_avoiding_pits(node_order_upstream, flow_receiver, node_z,
                        erosion_increment):
    n_nodes = node_order_upstream.size
    code="""
    int current_node;
    double elev_this_node_before;
    double elev_this_node_after;
    double elev_dstr_node_after;

    for (int i = 0; i < n_nodes; i++) {
        current_node = node_order_upstream[i];
        elev_this_node_before = node_z[current_node];
        elev_this_node_after = elev_this_node_before - erosion_increment[current_node];
        elev_dstr_node_after = elev_dstr[current_node] - erosion_increment[flow_receiver[current_node]];
        if (elev_this_node_after<elev_dstr_node_after) {
            erosion_increment[current_node] = (elev_this_node_before-elev_dstr_node_after)*0.999999;
        }
    }
    """
    weave.inline(code, ['n_nodes', 'node_order_upstream', 'node_z',
                        'erosion_increment', 'elev_dstr', 'flow_receiver'])


def erode_with_alpha(upstream_order_IDs, flow_receivers, alpha, z):
    code = """
    int current_node;
    int j;
    for (int i = 0; i < n_nodes; i++) {
        current_node = upstream_order_IDs[i];
        j = flow_receivers[current_node];
        if (current_node != j) {
            z[current_node] = (z[current_node] + alpha[current_node]*z[j])/(1.0+alpha[current_node]);
        }
    }
    """
    weave.inline(code, ['n_nodes', 'upstream_order_IDs', 'flow_receivers',
                        'z', 'alpha'])


def erode_with_link_alpha(upstream_order_IDs, flow_receivers,
                          alpha_by_flow_link_lengthtothenless1, n, z):
    code = """
    int current_node;
    int j;
    double current_z;
    double previous_z;
    double elev_diff;
    double elev_diff_tothenless1;
    for (int i = 0; i < n_nodes; i++) {
        current_node = upstream_order_IDs[i];
        j = flow_receivers[current_node];
        previous_z = z[current_node];
        if (current_node != j) {
            while (1) {
                elev_diff = previous_z-z[j];

                //this isn't defined if in some iterations the elev_diff
                // goes negative
                elev_diff_tothenless1 = pow(elev_diff, n-1.);
                current_z = previous_z - (previous_z - z[current_node] + alpha_by_flow_link_lengthtothenless1[current_node]*elev_diff_tothenless1*elev_diff)/(1.+n*alpha_by_flow_link_lengthtothenless1[current_node]*elev_diff_tothenless1);
                if (fabs((current_z - previous_z)/current_z) < 1.48e-08) {
                    break;
                }
                previous_z = current_z;
            }
            z[current_node] = current_z;
        }
    }
    """
    weave.inline(code, ['n_nodes', 'upstream_order_IDs', 'flow_receivers',
                        'z', 'alpha_by_flow_link_lengthtothenless1', 'n'],
                 headers=["<math.h>"])
