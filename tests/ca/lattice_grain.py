"""
Provide rules and states for a lattice-grain CTS model.

Created on Tue Feb  2 07:37:59 2016

@author: gtucker
"""

from landlab.ca.celllab_cts import Transition


def lattice_grain_node_states():
    """
    Create and return dict of states for lattice-grain model.

    Examples
    --------
    >>> nsd = lattice_grain_node_states()
    >>> len(nsd)
    9
    """
    ns_dict = {
        0: "empty",
        1: "moving up",
        2: "moving right and up",
        3: "moving right and down",
        4: "moving down",
        5: "moving left and down",
        6: "moving left and up",
        7: "rest",
        8: "wall",
    }
    return ns_dict


def lattice_grain_transition_list(g=0.0, f=0.0, motion=1.0, swap=False, callback=None):
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for simple granular mechanics model.

    Parameters
    ----------
    g : float (optional)
        Gravitational parameter. 1/g is the time taken for an initially
        stationary particle to fall one cell width. (Dimensions of 1/time)
    f : float (optional)
        Frictional parameter. Dimensions of 1/time. Probability per unit time
        of a frictional collision, in which one or both colliding particles
        come to a halt.
    motion : float (optional)
        Rate of motion (cells per time unit)
    swap : bool (optional)
        Whether to swap properties for motion transition; set to True if you
        want to track grain position via the 'propid' array.
    callback : function (optional)
        Function to be called after each transition. Used to implement custom
        handling of properties that evolve over time at each cell.

    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state transitions.

    Notes
    -----
    The transitions for this version of lattice gas have 11 pair-transition
    rules. The shorthand for the states is as follows:

        AR = air/empty
        IN = incoming particle (moving toward its neighbor)
        OU = outgoing particle (moving away from its neighbor)
        IO = incoming at an oblique angle
        OO = outgoing at an oblique angle
        RE = rest particle
        WA = wall particle
        op = oblique pair moving in opposite perpendicular direction
        sm = oblique pair moving in same perpendicular direction

    The 11 pairs with transitions are:

        1. AR-IN => IN-AR (move to empty particle)
        2. IN-IN => OO-OO-op (1/3 each dir), OU-OU (1/3) (head-on collision)
        3. IN-IO => OO-OU (oblique collision)
        4. IN-OO => IO-OU (oblique collision from behind)
        5. IN-OU => IO-OO (1/4 each of 2 directions) (collision from behind)
        6. IN-RE => RE-OU (1/3) RE-OO (1/3 each dir) (collision with rest)
        7. IN-WA => OU-WA (1/3) OO-WA (1/3 each dir) (wall collision)
        8. IO-IO-op => OO-OO-op (1/2 each dir) (glacing collision)
        9. IO-IO-sm => OO-OO-sm (30-degree collision)
        10. IO-RE => RE-OU (oblique collision with rest particle)
        11. IO-WA => OO-WA (oblique collision with wall)

    Examples
    --------
    >>> xnl = lattice_grain_transition_list()
    >>> len(xnl)  # 6 motion plus 117 elastic
    123
    >>> def fn():
    ...     return "test"
    ...
    >>> xnl = lattice_grain_transition_list(f=1.0, callback=fn)
    >>> len(xnl)  # 6 motion plus 87 frictional
    93
    >>> xnl[0].prop_update_fn()
    'test'
    >>> xnl = lattice_grain_transition_list(g=1.0, f=1.0, swap=True)
    >>> len(xnl)  # 6 motion plus 87 frictional plus 48 gravitational
    141
    >>> xnl[5].swap_properties
    True
    >>> xnl[6].swap_properties
    False
    """
    xn_list = []

    # Set elastic-response rate, and scale frictional-response rate to motion
    # rate
    p_elast = motion * (1.0 - f)  # rate of elastic (non-dissipative) collision
    f *= motion

    # Rule 1: Transitions for particle movement into an empty cell
    xn_list.append(Transition((1, 0, 0), (0, 1, 0), motion, "motion", swap, callback))
    xn_list.append(Transition((2, 0, 1), (0, 2, 1), motion, "motion", swap, callback))
    xn_list.append(Transition((3, 0, 2), (0, 3, 2), motion, "motion", swap, callback))
    xn_list.append(Transition((0, 4, 0), (4, 0, 0), motion, "motion", swap, callback))
    xn_list.append(Transition((0, 5, 1), (5, 0, 1), motion, "motion", swap, callback))
    xn_list.append(Transition((0, 6, 2), (6, 0, 2), motion, "motion", swap, callback))

    # Rule 2: Transitions for head-on collision: elastic
    if p_elast > 0.0:
        xn_list.append(
            Transition(
                (1, 4, 0), (4, 1, 0), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 4, 0), (3, 6, 0), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 4, 0), (5, 2, 0), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 5, 1), (5, 2, 1), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 5, 1), (4, 1, 1), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 5, 1), (6, 3, 1), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 6, 2), (6, 3, 2), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 6, 2), (1, 4, 2), p_elast / 3, "head-on collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 6, 2), (5, 2, 2), p_elast / 3, "head-on collision", False, callback
            )
        )

    # Rule 2: Transitions for head-on collision: frictional dissipation
    if f > 0.0:
        xn_list.append(
            Transition((1, 4, 0), (7, 7, 0), f, "head-on collision", False, callback)
        )
        xn_list.append(
            Transition((2, 5, 1), (7, 7, 1), f, "head-on collision", False, callback)
        )
        xn_list.append(
            Transition((3, 6, 2), (7, 7, 2), f, "head-on collision", False, callback)
        )

    # Rule 3: Transitions for oblique collision: elastic
    if p_elast > 0.0:
        xn_list.append(
            Transition(
                (1, 3, 0), (3, 1, 0), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 5, 0), (5, 1, 0), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 4, 0), (4, 2, 0), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (6, 4, 0), (4, 6, 0), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 4, 1), (4, 2, 1), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 6, 1), (6, 2, 1), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 5, 1), (5, 1, 1), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 5, 1), (5, 3, 1), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 1, 2), (1, 3, 2), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 5, 2), (5, 3, 2), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 6, 2), (6, 2, 2), p_elast, "oblique collision", False, callback
            )
        )
        xn_list.append(
            Transition(
                (4, 6, 2), (6, 4, 2), p_elast, "oblique collision", False, callback
            )
        )

    # Rule 3 frictional
    if f > 0.0:
        xn_list.append(
            Transition((1, 3, 0), (7, 7, 0), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((1, 5, 0), (7, 7, 0), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((2, 4, 0), (7, 7, 0), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((6, 4, 0), (7, 7, 0), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((2, 4, 1), (7, 7, 1), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((2, 6, 1), (7, 7, 1), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((1, 5, 1), (7, 7, 1), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((3, 5, 1), (7, 7, 1), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((3, 1, 2), (7, 7, 2), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((3, 5, 2), (7, 7, 2), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((2, 6, 2), (7, 7, 2), f, "oblique collision", False, callback)
        )
        xn_list.append(
            Transition((4, 6, 2), (7, 7, 2), f, "oblique collision", False, callback)
        )

    # Rule 4: Transitions for oblique-from-behind collisions
    if p_elast > 0.0:
        xn_list.append(
            Transition((1, 2, 0), (2, 1, 0), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((1, 6, 0), (6, 1, 0), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((3, 4, 0), (4, 3, 0), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((5, 4, 0), (4, 5, 0), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((2, 1, 1), (1, 2, 1), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((2, 3, 1), (3, 2, 1), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((4, 5, 1), (5, 4, 1), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((6, 5, 1), (5, 6, 1), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((3, 2, 2), (2, 3, 2), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((3, 4, 2), (4, 3, 2), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((1, 6, 2), (6, 1, 2), p_elast, "oblique", False, callback)
        )
        xn_list.append(
            Transition((5, 6, 2), (6, 5, 2), p_elast, "oblique", False, callback)
        )

    # Rule 4 frictional
    if f > 0.0:
        xn_list.append(Transition((1, 2, 0), (7, 1, 0), f, "oblique", False, callback))
        xn_list.append(Transition((1, 6, 0), (7, 1, 0), f, "oblique", False, callback))
        xn_list.append(Transition((3, 4, 0), (4, 7, 0), f, "oblique", False, callback))
        xn_list.append(Transition((5, 4, 0), (4, 7, 0), f, "oblique", False, callback))
        xn_list.append(Transition((2, 1, 1), (7, 2, 1), f, "oblique", False, callback))
        xn_list.append(Transition((2, 3, 1), (7, 2, 1), f, "oblique", False, callback))
        xn_list.append(Transition((4, 5, 1), (5, 7, 1), f, "oblique", False, callback))
        xn_list.append(Transition((6, 5, 1), (5, 7, 1), f, "oblique", False, callback))
        xn_list.append(Transition((3, 2, 2), (7, 3, 2), f, "oblique", False, callback))
        xn_list.append(Transition((3, 4, 2), (7, 3, 2), f, "oblique", False, callback))
        xn_list.append(Transition((1, 6, 2), (6, 7, 2), f, "oblique", False, callback))
        xn_list.append(Transition((5, 6, 2), (6, 7, 2), f, "oblique", False, callback))

    # Rule 5: Transitions for direct-from-behind collisions
    if p_elast > 0.0:
        xn_list.append(
            Transition((1, 1, 0), (2, 6, 0), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((1, 1, 0), (6, 2, 0), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((4, 4, 0), (3, 5, 0), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((4, 4, 0), (5, 3, 0), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((2, 2, 1), (1, 3, 1), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((2, 2, 1), (3, 1, 1), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((5, 5, 1), (4, 6, 1), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((5, 5, 1), (6, 4, 1), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((3, 3, 2), (2, 4, 2), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((3, 3, 2), (4, 2, 2), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((6, 6, 2), (1, 5, 2), p_elast / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((6, 6, 2), (5, 1, 2), p_elast / 4, "behind", False, callback)
        )

    # Rule 5 frictional
    if f > 0.0:
        xn_list.append(
            Transition((1, 1, 0), (7, 1, 0), f / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((4, 4, 0), (4, 7, 0), f / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((2, 2, 1), (7, 2, 1), f / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((5, 5, 1), (5, 7, 1), f / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((3, 3, 2), (7, 3, 2), f / 4, "behind", False, callback)
        )
        xn_list.append(
            Transition((6, 6, 2), (6, 7, 2), f / 4, "behind", False, callback)
        )

    # Rule 6: Transitions for direct collision with stationary (resting) particle
    if p_elast > 0.0:
        xn_list.append(
            Transition((1, 7, 0), (7, 1, 0), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((1, 7, 0), (7, 2, 0), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((1, 7, 0), (7, 6, 0), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 4, 0), (4, 7, 0), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 4, 0), (3, 7, 0), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 4, 0), (5, 7, 0), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((2, 7, 1), (7, 2, 1), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((2, 7, 1), (7, 1, 1), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((2, 7, 1), (7, 3, 1), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 5, 1), (5, 7, 1), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 5, 1), (4, 7, 1), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 5, 1), (6, 7, 1), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((3, 7, 2), (7, 3, 2), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((3, 7, 2), (7, 2, 2), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((3, 7, 2), (7, 4, 2), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 6, 2), (6, 7, 2), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 6, 2), (1, 7, 2), p_elast / 3.0, "rest", False, callback)
        )
        xn_list.append(
            Transition((7, 6, 2), (5, 7, 2), p_elast / 3.0, "rest", False, callback)
        )

    # Rule 6 frictionl
    if f > 0.0:
        xn_list.append(Transition((1, 7, 0), (7, 7, 0), f, "rest", False, callback))
        xn_list.append(Transition((7, 4, 0), (7, 7, 0), f, "rest", False, callback))
        xn_list.append(Transition((2, 7, 1), (7, 7, 1), f, "rest", False, callback))
        xn_list.append(Transition((7, 5, 1), (7, 7, 1), f, "rest", False, callback))
        xn_list.append(Transition((3, 7, 2), (7, 7, 2), f, "rest", False, callback))
        xn_list.append(Transition((7, 6, 2), (7, 7, 2), f, "rest", False, callback))

    # Rule 7: Transitions for wall impact
    if p_elast > 0.0:
        xn_list.append(
            Transition(
                (1, 8, 0), (4, 8, 0), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 8, 0), (3, 8, 0), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 8, 0), (5, 8, 0), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 8, 1), (5, 8, 1), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 8, 1), (4, 8, 1), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 8, 1), (6, 8, 1), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 8, 2), (6, 8, 2), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 8, 2), (5, 8, 2), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 8, 2), (1, 8, 2), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 4, 0), (8, 1, 0), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 4, 0), (8, 6, 0), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 4, 0), (8, 2, 0), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 5, 1), (8, 1, 1), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 5, 1), (8, 2, 1), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 5, 1), (8, 3, 1), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 6, 2), (8, 2, 2), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 6, 2), (8, 3, 2), p_elast / 3, "wall rebound", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 6, 2), (8, 4, 2), p_elast / 3, "wall rebound", False, callback
            )
        )

    # Rule 7 frictional
    if f > 0.0:
        xn_list.append(
            Transition((1, 8, 0), (7, 8, 0), f, "wall rebound", False, callback)
        )
        xn_list.append(
            Transition((2, 8, 1), (7, 8, 1), f, "wall rebound", False, callback)
        )
        xn_list.append(
            Transition((3, 8, 2), (7, 8, 2), f, "wall rebound", False, callback)
        )
        xn_list.append(
            Transition((8, 4, 0), (8, 7, 0), f, "wall rebound", False, callback)
        )
        xn_list.append(
            Transition((8, 5, 1), (8, 7, 1), f, "wall rebound", False, callback)
        )
        xn_list.append(
            Transition((8, 6, 2), (8, 7, 2), f, "wall rebound", False, callback)
        )

    # Rule 8: Transitions for glancing oblique collision
    if p_elast > 0.0:
        xn_list.append(
            Transition((2, 5, 0), (3, 6, 0), p_elast, "glancing", False, callback)
        )
        xn_list.append(
            Transition((6, 3, 0), (5, 2, 0), p_elast, "glancing", False, callback)
        )
        xn_list.append(
            Transition((3, 6, 1), (4, 1, 1), p_elast, "glancing", False, callback)
        )
        xn_list.append(
            Transition((1, 4, 1), (6, 3, 1), p_elast, "glancing", False, callback)
        )
        xn_list.append(
            Transition((4, 1, 2), (5, 2, 2), p_elast, "glancing", False, callback)
        )
        xn_list.append(
            Transition((2, 5, 2), (1, 4, 2), p_elast, "glancing", False, callback)
        )

    # Rule 8 frictional
    if f > 0.0:
        xn_list.append(Transition((2, 5, 0), (7, 7, 0), f, "glancing", False, callback))
        xn_list.append(Transition((6, 3, 0), (7, 7, 0), f, "glancing", False, callback))
        xn_list.append(Transition((3, 6, 1), (7, 7, 1), f, "glancing", False, callback))
        xn_list.append(Transition((1, 4, 1), (7, 7, 1), f, "glancing", False, callback))
        xn_list.append(Transition((4, 1, 2), (7, 7, 2), f, "glancing", False, callback))
        xn_list.append(Transition((2, 5, 2), (7, 7, 2), f, "glancing", False, callback))

    # Rule 9: Transitions for "near-on" collisions
    if p_elast > 0.0:
        xn_list.append(
            Transition((6, 5, 0), (5, 6, 0), p_elast, "near-on", False, callback)
        )
        xn_list.append(
            Transition((2, 3, 0), (3, 2, 0), p_elast, "near-on", False, callback)
        )
        xn_list.append(
            Transition((1, 6, 1), (6, 1, 1), p_elast, "near-on", False, callback)
        )
        xn_list.append(
            Transition((3, 4, 1), (4, 3, 1), p_elast, "near-on", False, callback)
        )
        xn_list.append(
            Transition((2, 1, 2), (1, 2, 2), p_elast, "near-on", False, callback)
        )
        xn_list.append(
            Transition((4, 5, 2), (5, 4, 2), p_elast, "near-on", False, callback)
        )

    # Rule 9 frictional
    if f > 0.0:
        xn_list.append(
            Transition((6, 5, 0), (7, 6, 0), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((6, 5, 0), (5, 7, 0), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((2, 3, 0), (7, 2, 0), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((2, 3, 0), (3, 7, 0), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((1, 6, 1), (7, 1, 1), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((1, 6, 1), (6, 7, 1), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((3, 4, 1), (7, 3, 1), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((3, 4, 1), (4, 7, 1), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((2, 1, 2), (7, 2, 2), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((2, 1, 2), (1, 7, 2), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((4, 5, 2), (7, 4, 2), f / 2, "near-on", False, callback)
        )
        xn_list.append(
            Transition((4, 5, 2), (5, 7, 2), f / 2, "near-on", False, callback)
        )

    # Rule 10: Transitions for oblique collision with rest particle
    if p_elast > 0.0:
        xn_list.append(
            Transition(
                (2, 7, 0), (7, 1, 0), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (6, 7, 0), (7, 1, 0), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (7, 3, 0), (4, 7, 0), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (7, 5, 0), (4, 7, 0), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 7, 1), (7, 2, 1), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 7, 1), (7, 2, 1), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (7, 6, 1), (5, 7, 1), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (7, 4, 1), (5, 7, 1), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (4, 7, 2), (7, 3, 2), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 7, 2), (7, 3, 2), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (7, 5, 2), (6, 7, 2), p_elast, "oblique with rest", False, callback
            )
        )
        xn_list.append(
            Transition(
                (7, 1, 2), (6, 7, 2), p_elast, "oblique with rest", False, callback
            )
        )

    # Rule 10 frictional
    if f > 0.0:
        xn_list.append(
            Transition((2, 7, 0), (7, 7, 0), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((6, 7, 0), (7, 7, 0), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((7, 3, 0), (7, 7, 0), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((7, 5, 0), (7, 7, 0), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((3, 7, 1), (7, 7, 1), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((1, 7, 1), (7, 7, 1), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((7, 6, 1), (7, 7, 1), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((7, 4, 1), (7, 7, 1), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((4, 7, 2), (7, 7, 2), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((2, 7, 2), (7, 7, 2), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((7, 5, 2), (7, 7, 2), f, "oblique with rest", False, callback)
        )
        xn_list.append(
            Transition((7, 1, 2), (7, 7, 2), f, "oblique with rest", False, callback)
        )

    # Rule 11: Transitions for oblique collision with wall particle
    if p_elast > 0.0:
        xn_list.append(
            Transition(
                (2, 8, 0), (3, 8, 0), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (6, 8, 0), (5, 8, 0), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 3, 0), (8, 2, 0), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 5, 0), (8, 6, 0), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (1, 8, 1), (6, 8, 1), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (3, 8, 1), (4, 8, 1), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 4, 1), (8, 3, 1), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 6, 1), (8, 1, 1), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (4, 8, 2), (5, 8, 2), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (2, 8, 2), (1, 8, 2), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 1, 2), (8, 2, 2), p_elast, "oblique with wall", False, callback
            )
        )
        xn_list.append(
            Transition(
                (8, 5, 2), (8, 4, 2), p_elast, "oblique with wall", False, callback
            )
        )

    # Rule 11 frictional
    if f > 0.0:
        xn_list.append(
            Transition((2, 8, 0), (7, 8, 0), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((6, 8, 0), (7, 8, 0), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((8, 3, 0), (8, 7, 0), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((8, 5, 0), (8, 7, 0), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((1, 8, 1), (7, 8, 1), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((3, 8, 1), (7, 8, 1), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((8, 4, 1), (8, 7, 1), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((8, 6, 1), (8, 7, 1), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((4, 8, 2), (7, 8, 2), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((2, 8, 2), (7, 8, 2), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((8, 1, 2), (8, 7, 2), f, "oblique with wall", False, callback)
        )
        xn_list.append(
            Transition((8, 5, 2), (8, 7, 2), f, "oblique with wall", False, callback)
        )

    # Gravity rule 1: rising particles become rest particles
    if g > 0.0:
        xn_list.append(
            Transition((0, 1, 0), (0, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((1, 1, 0), (1, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((2, 1, 0), (2, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((3, 1, 0), (3, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((4, 1, 0), (4, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((5, 1, 0), (5, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((6, 1, 0), (6, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((7, 1, 0), (7, 7, 0), g, "gravity 1", False, callback)
        )
        xn_list.append(
            Transition((8, 1, 0), (8, 7, 0), g, "gravity 1", False, callback)
        )

    # Gravity rule 2: resting particles become falling particles (if above air)
    if g > 0.0:
        xn_list.append(
            Transition((0, 7, 0), (0, 4, 0), g, "gravity 2", False, callback)
        )

    # Gravity rule 3: up/sideways particles become down/sideways particles
    if g > 0.0:
        xn_list.append(
            Transition((0, 2, 0), (0, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((1, 2, 0), (1, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((2, 2, 0), (2, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((3, 2, 0), (3, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((4, 2, 0), (4, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((5, 2, 0), (5, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((6, 2, 0), (6, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((7, 2, 0), (7, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((8, 2, 0), (8, 3, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((0, 6, 0), (0, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((1, 6, 0), (1, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((2, 6, 0), (2, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((3, 6, 0), (3, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((4, 6, 0), (4, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((5, 6, 0), (5, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((6, 6, 0), (6, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((7, 6, 0), (7, 5, 0), g, "gravity 3", False, callback)
        )
        xn_list.append(
            Transition((8, 6, 0), (8, 5, 0), g, "gravity 3", False, callback)
        )

    # Gravity rule 4: down/side to straight down
    if g > 0.0:
        xn_list.append(
            Transition((0, 3, 0), (0, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((1, 3, 0), (1, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((2, 3, 0), (2, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((3, 3, 0), (3, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((4, 3, 0), (4, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((5, 3, 0), (5, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((6, 3, 0), (6, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((7, 3, 0), (7, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((8, 3, 0), (8, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((0, 5, 0), (0, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((1, 5, 0), (1, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((2, 5, 0), (2, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((3, 5, 0), (3, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((4, 5, 0), (4, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((5, 5, 0), (5, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((6, 5, 0), (6, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((7, 5, 0), (7, 4, 0), g, "gravity 4", False, callback)
        )
        xn_list.append(
            Transition((8, 5, 0), (8, 4, 0), g, "gravity 4", False, callback)
        )

    # Gravity rule for lateral destabilization (represents grain
    # motion above angle of repose on sloping surface)
    if g > 0.0:
        xn_list.append(
            Transition((7, 0, 2), (3, 0, 2), g / 2.0, "gravity", False, callback)
        )
        xn_list.append(
            Transition((0, 7, 1), (0, 5, 1), g / 2.0, "gravity", False, callback)
        )

    return xn_list
