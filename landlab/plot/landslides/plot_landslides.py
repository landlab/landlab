# -*- coding: utf-8 -*-
"""plot_landslides.py
Created on Mon May 9 2020
@author: christinab
"""
import matplotlib.pyplot as plt
from landlab.plot import imshow_grid
import numpy as np

def test_field_initialized_to_range(ls_prob,input_var_names,unit_default_value):        
    """Testing if the output fields are initialized within range of default.
    """
    for name in input_var_names:
        
        print("Input: (Min, Max) range of {value} model inputs: {pars}".format(value=name,pars=unit_default_value[name]))
        print("Ouput: {value} default value initialized for each node:".format(value=name))
        field = ls_prob.grid["node"][name]
        print(field)
        
        field_min = ls_prob.grid["node"][name].min()
        field_max = ls_prob.grid["node"][name].max()
        unit_min = unit_default_value[name][0]
        unit_max = unit_default_value[name][1]

        assert field_min >= unit_min
        assert field_max <= unit_max 

        print("Test: {value} within range. Passes test.".format(value=name))
        print("")
        
def print_nodevalues(ls_prob):
    """Print out the input range and resulting node values for 
    all landslide component values on the node.
    ----------
    Input: Landlab landslide model instance  """
    
    for name in ls_prob.grid["node"]:
        field_min = ls_prob.grid["node"][name].min()
        field_max = ls_prob.grid["node"][name].max()

        print_minmax = print("{value} (Min, Max) range of ({value1},{value2}) values on node:".format(value=name,value1=field_min,value2=field_max))
        print_value = print("All {value} default value initialized for each node:".format(value=name))
        field = ls_prob.grid["node"][name]
        print(field)

        print("")

def print_list_nodevalues(ls_prob,values_list):
    """Print out the input range and resulting node values for 
    all landslide component values on the node.
    ----------
    Input: Landlab landslide model instance  """
    
    for name in values_list:
        field_min = ls_prob.grid["node"][name].min()
        field_max = ls_prob.grid["node"][name].max()

        print_minmax = print("{value} (Min, Max) range of ({value1},{value2}) values on node:".format(value=name,value1=field_min,value2=field_max))
        print_value = print("All {value} default value initialized for each node:".format(value=name))
        field = ls_prob.grid["node"][name]
        print(field)

        print("")

def plot_landslide_4variables(littleplots,
                              inputgrid1,inputgrid2,inputgrid3,inputgrid4,
                              variable,maxval,plot_subtitle,cmap,scalelabel):
    """Return figure with four subplots given landslide variables.
    Parameters
    ----------
    name : littleplots
    List of subplot matrix locations: row, column, index
    Examples: 
    littleplots = [221,222,223,224]
    littleplots = [411,412,413,414]
    littleplots = [141,142,143,144]
    """
    
    fig = plt.figure('Plot Component Inputs/Outputs')

    ax1=[]
    ax2=[]
    ax3=[]
    ax4=[]
    for eachgrid, plotvar, maxval,plot_name,cmap, cbar_name, littleplots,axeslist in zip([inputgrid1,inputgrid2,inputgrid3,inputgrid4],
                                                      [inputgrid1.at_node[variable[0]],
                                                       inputgrid2.at_node[variable[1]],
                                                       inputgrid3.at_node[variable[2]],
                                                       inputgrid4.at_node[variable[3]]],
                                                       maxval,        
                                                       plot_subtitle,
                                                       [cmap[0],cmap[1],cmap[2],cmap[3]],
                                                       [scalelabel[0],scalelabel[1],scalelabel[2],scalelabel[3]],
                                                       littleplots,
                                                       [ax1,ax2,ax3,ax4]):


            axeslist = fig.add_subplot(littleplots)
            axeslist.xaxis.set_visible(False)

            imshow_grid(eachgrid, plotvar, plot_name=plot_name, cmap=cmap,
                     grid_units=('coordinates', 'coordinates'), shrink=0.75,limits=(0, maxval),
                     var_name=cbar_name)

    fig.set_size_inches(10,10)
    return None

    
def scenario_unit_explorer(rw,a,T,theta,hs):
    print("Given relative wetness = {value}".format(value=rw))
    #print("Mean topographic__specific_contributing_area= {value_mean} ".format(value_mean=a.mean()))
    
    Recharge = ((rw * (T * theta )) / a)*1000 #mm/day

    Remin_value = round(Recharge.min())
    Remean = round(Recharge.mean())
    Restandard_deviation =  round(Recharge.std())
    Remax_value =  round(Recharge.max())
    rw_r = Recharge * a / (T * theta)

    Scenario_R=[Remin_value,Remax_value,Remean,Restandard_deviation]

    
    De_sat_threshold = 0.001
    Depth = hs - rw * hs
    rw_d =  (hs - Depth) / (hs-De_sat_threshold)
    hw = hs - Depth
 
    Demin_value = round(Depth.min(),2) #assumes unit test with soil thickness 1 m
    Demean = round(Depth.mean(),2)
    Destandard_deviation =  round(Depth.std(),2)
    Demax_value =  round(Depth.max(),2)

    Scenario_D=[Demin_value.min(),Demax_value.max(),Demean.mean(),Destandard_deviation.mean()]
  
    #print("Recharge min, max, mean, std  {value1} [min,max,mean,std] meters/day".format(value1=Scenario_R))
    #print("Depth min, max, mean, std  {value1} [min,max,mean,std] meters".format(value1=Scenario_D))
    print("Unit Recharge Parameters mm/day:")
    print("Remin = {value} ".format(value=Remin_value))
    print("Remean = {value} ".format(value=Remean))
    print("Restandard_deviation = {value} ".format(value=Restandard_deviation))
    print("Remax = {value} ".format(value=Remax_value))
    #print("Relative Wetness from Recharge= {value} ".format(value=rw_r))
    print("")
    #print("Depth to groundwater= {value} ".format(value=np.round(Depth,2)))
    #print("Height of water= {value} ".format(value=np.round(hw,2)))
    #print("Relative Wetness from Depth= {value} ".format(value=np.round(rw_d,2)))
    print("")
    print("Unit Recharge Parameters mm/day:")
    print("Demin_value = {value} ".format(value=Scenario_D[0]))
    print("Demax_value = {value} ".format(value=Scenario_D[1]))
    print("Demean = {value} ".format(value=Scenario_D[2]))
    print("Destandard_deviation = {value} ".format(value=Scenario_D[3]))

        
    return Scenario_R, Scenario_D   
 