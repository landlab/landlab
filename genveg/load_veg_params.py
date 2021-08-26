"""
Load vegetation parameters from formatted Excel file for use
in GenVeg, defines VegParams class
"""

import pandas as pd

#Create veg_params class type for input data and calculated params
class VegParams:
    "This is a basic vegetation type"
    def __init__(self,input):
        ##########################
        #Vegetation community info
        ##########################
        self.name=input[0]
        self.type=input[1]
        #Photosynthesis variables
        self.ph_type=input[3]
        #Growing season info
        self.gs_start=input[4]
        self.gs_end=input[5]
        self.senes=input[6]
        #Respiration coefficients
        self.res_co=input[58]
        #Glucose requirements
        self.glu_req=input[59]
        #Light extinction coefficient
        self.le_k=input[9]
        #Light half-saturation constant
        self.hi=input[10]
        #Maximal gross photosynthesis
        self.p_max=input[11]
        #######################
        #Plant size constraints
        #######################
        #Max plants per unit area
        self.pl_dens_max=input[13]
        #Max stems per plant
        self.stems_max=input[14]
        #Max stem height
        self.stem_height=input[15]
        #Max stem mass
        self.stem_mass=input[16]
        #Allocation of biomass between roots, stems, leaves
        self.allocate=input[60]
        #Stem cross-sectional area
        self.tca=input[18]
        #Minimum leaf mass
        self.leaf_min=input[19]
        #Minimum root mass
        self.root_min=input[20]
        #####################
        #Dispersal parameters
        #####################
        #Maximum dispersal distance per growing season
        self.disp_dist=input[22]
        #Minimum plant size for dispersal
        self.disp_size=input[23]
        #Dispersal carbohydrate cost
        self.disp_cost=input[24]
        ########################
        #Colonization parameters
        ########################
        self.col_prob=input[26]
        self.col_n=input[27]
        ##########################
        #Winter storage parameters
        ##########################
        self.wint_die=input[29]
        self.wint_store=input[30]

    #Calculate root:shoot ratio    
    def rsratio(self):
        return self.allocate[0]/(self.allocate[1]+self.allocate[2])

    #Calculate the mass of a stem per unit length
    def stem_mass_unit(self):
        return self.stem_mass/self.stem_height
    
    #Calculate the maximum stem density per unit area
    def stem_dens_max(self):
        return self.pl_dens_max*self.stems_max

    #Calculate the maximum aboveground biomass per unit area
    def ag_mass_max(self):
        return self.stem_dens_max()*self.stem_mass*(1+(self.allocate[2]/self.allocate[1]))
    
    #Calculate the maximum belowground biomass per unit area
    def bg_mass_max(self):
        return self.ag_mass_max()*self.rsratio()
    
    #Calculate growing season length
    def gs_length(self):
        return self.gs_end-self.gs_start+1

def load_veg_params(file):

    #file='GenVeg_params_inputs.xlsx'
    xlin=pd.ExcelFile(file)
    coms=xlin.sheet_names
    for i in coms:
        df_in=xlin.parse(i, usecols='C:H')
        samp_input=df_in['Values'].tolist()

        samp_input.append((df_in.Values[7],df_in.Val1[7],df_in.Val2[7]))
        samp_input.append((df_in.Values[8],df_in.Val1[8],df_in.Val2[8]))
        samp_input.append((df_in.Values[17],df_in.Val1[17],df_in.Val2[17]))
        x[i]=VegParams(samp_input)
    return x
