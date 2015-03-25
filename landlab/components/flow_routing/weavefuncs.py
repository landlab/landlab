from scipy import weave


def adjust_flow_receivers(fromnode, tonode, elev, link_slope, active_links,
                          receiver, receiver_link, steepest_slope):
    n_nodes = fromnode.size

    code="""
    int f;
    int t;
    for (int i=0; i < n_nodes; i++) {
        f = fromnode[i];
        t = tonode[i];
        if (elev[f] > elev[t] && link_slope[i] > steepest_slope[f]) {
            receiver[f] = t;
            steepest_slope[f] = link_slope[i];
            receiver_link[f] = active_links[i];
        }
        else if (elev[t] > elev[f] && - link_slope[i] > steepest_slope[t]) {
            receiver[t] = f;
            steepest_slope[t] = -link_slope[i];
            receiver_link[t] = active_links[i];
        }
    }
    """

    weave.inline(code, ['n_nodes', 'fromnode', 'tonode', 'elev',
                        'link_slope', 'steepest_slope', 'receiver',
                        'receiver_link', 'active_links'])
