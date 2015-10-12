#!/usr/env/python

"""
delta_liang_ca.py

Modifies Man Liang's reduced complexity delta evolution model, implementing a
pair-based scheme after Narteau.
Possible issues using a D4 not D8 scheme. We'll see when we run.

DEJH, Aug 2014
"""

import time
import numpy
from landlab import RasterModelGrid
from landlab.components.linkca.link_ca import LinkCellularAutomaton, Transition, CAPlotter

def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for a weathering model.

    We employ a scheme where sediment *can* move independently of water
    precipitons, but once it is caught up by a precipiton, it tends to stay
    with it until ejected to a floodplain, dropped by the flow, or the
    precipiton reaches the sea.

    The states and transitions are as follows (note: X-X means "horizontal"
    pair, X/X means "vertical" pair with first item beneath second):

    Node states
    -----------
    0 ocean
    1 channel
    2 floodplain
    3 channel with water packet
    4 channel with sediment packet
    5 channel with water and sediment packet
    6 floodplain with sediment packet
    7 abandoned channel

    Pair state      Transition to       Process     Rate
    ----------      -------------       -------     ----
    (0-0)
    (0-1)
    (0-2)
    (0-3)           (1-1)               Channel extends
                    (0-1)               Small chance channelization doesn't occur
    (0-4)           (2-1)               Deposition of channel sed as a floodplain bar
    (0-5)           (4-1)               Channel extends,uses up its precipiton, sed packet moves along with it
    (0-6)           (2-2)               Floodplain grows into ocean
    (0-7)           (0-0)               Abandoned channel inundated by ocean (could alter this rule later)
    (1-1)
    (1-2)
    (1-3)           (3-1)               Water passes along channel
    (1-4)           (4-1)               Sediment passes along channel
    (1-5)           (3-4)               Sediment deposited by precipiton
                    (5-1)               precipiton moves downstream carrying its sediment load
    (1-6)           (4-2)               Sediment enters channel
    (1-7)
    (2-2)
    (2-3)
    (2-4)           (6-1)               Sediment packet moves from channel to floodplain
    (2-5)           (6-3)               Sediment moves from wet channel to floodplain (flood deposition)
    (2-6)           (6-2)               Sediment moves across floodplain
    (2-7)
    (3-3)
    (3-4)           (1-5)               Water envelops sediment packet
                    (5-1)               Sediment packet migrates into precipiton
    (3-5)           (5-3)               Sed packet gets passed downstream in flow
    (3-6)           (5-2)               Sed cannibalized from floodplain into channel
    (3-7)
    (4-4)
    (4-5)           (5-4)               Flow moves downstream over a heavily sedimented bed
    (4-6)
    (4-7)           (1-2)               Plug forms (will be an abandoned channel behind it)
    (5-5)
    (5-6)
    (5-7)           (3-2)               Plug forms (will be an abandoned channel behind it)
    (6-6)
    (6-7)           (2-2)               Abandoned channel fills with floodplain sed
    (7-7)

    """
    xn_list = []

    subaerial_sed_mobility = 1.
    precipiton_mobility = 1.
    jump_to_plain_rate = subaerial_sed_mobility
    backwater_rate = precipiton_mobility #This is the (3-5)->(5-3) rate. set it lower to simulate backwater effects??


    xn_list.append( Transition((0,3,0), (1,0,0), ) ) #there won't be direction here... probably (Man does use a global direction wrighting, gamma)
