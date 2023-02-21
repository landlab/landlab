"""
Load veg veg_params
May need to change this to output a landlab configured input file
Want to retain formatted Excel file due to type of data required
"""
import numpy as np
import pandas as pd
import pathlib
import yaml
from scipy.optimize import curve_fit

#Create veg_params class type for input data and calculated params
class VegParams:
    """
    Load vegetation parameters from formatted Excel file
    into a Pandas dateframe that is converted to parameter 
    dictionaries stored in a yaml for landlab input

    Parameters
    ----------
    fpath: Pathfile object, file path where input file is located. Must be structured Excel or csv file or yaml.
           If blank, model will assume species is corn.
    processes: a list of vegetation processes to initialize parameters for 'plantsize','dispersal','mortality','colonization'. 
               If blank, model will run basic growth only.
    outfile: optional input to allow for custom file name
    """
    def __init__(
        self,
        fpath='None',
        outfile='veg_params.yml',
        processes=[],
        vegparams={}
    ):
        if fpath=='None':
            self.veg_params={'Corn': {
                    'plant_factors':{
                        'species':'Corn',
                        'growth_habit': 'forb_herb',
                        'monocot_dicot': 'monocot',
                        'angio_gymno': 'angiosperm',
                        'annual_perennial': 'annual',
                        'leaf_retention': 'deciduous',
                        'growth_form':'single_stem',
                        'shape':'erect',
                        'ptype':'C3'
                    },
                    'grow_params': {
                        'growing_season_start': 91,
                        'growing_season_end': 290,
                        'senescence_start': 228,
                        'respiration_coefficient': [0.015,0.015,0.03],
                        'glucose_requirement': [1.444,1.513,1.463],
                        'k_light_extinct':0.02,
                        'light_half_sat':9,
                        'p_max':0.055,
                        'root_to_leaf_coeffs': [0.031,0.951,0],
                        'root_to_stem_coeffs': [-0.107, 1.098, 0.0216],
                        'plant_part_min':[0.01,0.1,0.5],
                        'plant_part_max':[6,25,30]
                    }
                }
            }
            if 'plantsize' in processes:
                self.size_params={
                    'max_plant_density': 1,
                    'max_n_stems': 3,
                    'max_height_stem': 2.5,
                    'max_mass_stem': 72,
                    'total_cs_area_stems': 0.231                    
                }
                if 'dispersion' in processes:
                    self.disp_params={
                        'max_dist_dispersal': 2,
                        'disp_size_rat': 0.5,
                        'disp_cost': 0
                    }
                else: self.disp_params={}
                self.veg_params['Corn']['dispersal_params']={**self.dispersal_params}
            else: self.size_params={}
            self.veg_params['Corn']['size_params']={**self.size_params}
            if 'colonize' in processes:
                self.col_params={
                    'col_prob': 0.01,
                    'col_dt': 365
                }
            else: self.col_params={}
            self.veg_params['Corn']['col_params']={**self.col_params}
            if 'mortality' in processes:
                self.mort_params={
                    'mort_factor_1': 'Mortality factor',
                    'mort_factor_1_duration': 365,
                    'mort_factor_1_coeffs': [0,0]
                }
            else: self.mort_params={}
            self.vegparams['Corn']['mort_params']={**self.mort_params}
        else: 
            ispathvalid=fpath.is_file()   
            if ispathvalid==False:
                raise ValueError('File path is not valid')        
            self.fpath=fpath
        #add check for file extension
            exten=pathlib.Path(self.fpath).suffix
            if exten == 'yml':
                print('File already in correct file format. Use Landlab load_params function.')
                pass
            else:
                if 'xls' in exten:          
                #Read Excel file data into dataframe
                    xlin=pd.ExcelFile(fpath)
                    coms=xlin.sheet_names
                    x=pd.DataFrame()
                    #Create list of values for a community on each tab
                    for i in coms:    
                        df_in=xlin.parse(i, usecols='B,C,D')
                        df_in.set_index(['Variable Name','Descriptor'], inplace=True)
                        #Define plant factor keys and create plant factor dictionary
                        factor_keys=[
                            'species',
                            'growth_habit',
                            'monocot_dicot',
                            'angio_gymno',
                            'duration',
                            'leaf_retention',
                            'growth_form',
                            'shape',
                            'p_type'
                        ]
                        factor=self._makedict(df_in, factor_keys, 'plant_factors')    
                        
                        duration_keys=[
                            'growing_season_start',
                            'growing_season_end',
                            'senescence_start',
                        ]
                        duration=self._makedict(df_in, duration_keys,'duration_params')

                        #Define growth parameter keys and create growthparams dictionary                   
                        growth_keys=[
                            'respiration_coefficient',
                            'glucose_requirement',
                            'k_light_extinct',
                            'light_half_sat',
                            'p_max',
                            'plant_part_min',
                            'plant_part_max',
                            'root_to_leaf_coeffs',
                            'root_to_stem_coeffs'
                        ]

                        #Replace null values for coefficients with Poorter-derived coefficients if necessary
                        woody_herb=('herb','woody')[factor['plant_factors']['growth_habit']=='shrub']
                        opt_2=(factor['plant_factors']['monocot_dicot'], factor['plant_factors']['angio_gymno'])[woody_herb=='woody']
                        options={
                            'woody':{
                                'angiosperm':{'root_to_leaf_coeffs':[0.090,0.889,-0.0254],'root_to_stem_coeffs':[-0.097,1.071,0.0179]},
                                'gymnosperm':{'root_to_leaf_coeffs':[0.243,0.924,-0.0282],'root_to_stem_coeffs':[-0.070,1.236,-0.0186]}
                                },
                            'herb':{
                                'monocot':{'root_to_leaf_coeffs':[0.031,0.951,0],'root_to_stem_coeffs':[-0.107,1.098,0.0216]},
                                'dicot':{'root_to_leaf_coeffs':[0.259,0.916,0],'root_to_stem_coeffs':[-0.111,1.029,0]}
                                }
                            }
                        df_fill=options[woody_herb][opt_2]
                        if df_in.loc['root_to_leaf_coeffs'].isnull().values.any():
                            df_in.loc['root_to_leaf_coeffs','Values']=df_fill['root_to_leaf_coeffs']
                        if df_in.loc['root_to_stem_coeffs'].isnull().values.any():
                            df_in.loc['root_to_stem_coeffs','Values']=df_fill['root_to_stem_coeffs']
                        
                        grow=self._makedict(df_in, growth_keys, 'grow_params')
                        multi_part_keys=['glucose_requirement', 'respiration_coefficient','plant_part_min','plant_part_max']
                        for key in multi_part_keys:
                            grow['grow_params'][key]=self.replace_part_list_with_dict(grow['grow_params'][key])
                        grow['grow_params']['total_min_biomass']=sum(grow['grow_params']['plant_part_min'].values())
                        grow['grow_params']['total_max_biomass']=sum(grow['grow_params']['plant_part_max'].values())
                        grow['grow_params']['growth_min_biomass']=grow['grow_params']['total_min_biomass']-grow['grow_params']['plant_part_min']['storage']-grow['grow_params']['plant_part_min']['reproductive']
                        grow['grow_params']['growth_max_biomass']=grow['grow_params']['total_max_biomass']-grow['grow_params']['plant_part_max']['storage']-grow['grow_params']['plant_part_max']['reproductive']
                        #If plantsize is required process, define plant size parameter keys and create sizeparams dictionary
                        if 'plantsize' in processes:
                            size_keys=[
                                'max_plant_density',
                                'max_n_stems',
                                'max_height_stem',
                                'total_cs_area_stems'
                            ]
                            size=self._makedict(df_in, size_keys, 'size_params')
                            #If dispersion is required process, define dispersion parameter keys and create dispparams dictionary
                            if 'dispersal' in processes:
                                disp_keys=[
                                    'reproduction_start',
                                    'max_dist_dispersal',
                                    'min_size_dispersal',
                                    'carb_cost_dispersal'
                                ]
                                disp=self._makedict(df_in, disp_keys, 'dispersal_params')
                            else:
                                disp={}
                        else:
                            size={}
                        #If colonization is required process, define coloniation parameters keys and create colparams dictionary
                        if 'colonize' in processes:
                            col_keys=[
                                'prob_colonization',
                                'time_to_colonization'
                            ]
                            col=self._makedict(df_in, col_keys, 'col_params')
                        else:
                            col={}
                        #If mortality is required process, define mortality parameter keys and create mortparams dictionary
                        if 'mortality' in processes:
                            print(df_in)
                            mort_df=df_in.dropna()
                            mort_df_group=mort_df.groupby(['Descriptor','Variable Name']).agg(pd.Series.tolist)
                            factors=mort_df_group.index.levels[0]
                            mort_dict={}
                            for factor in factors:
                                mort_vars=mort_df_group.xs(factor)
                                names=mort_vars['Values'].loc['name'][0]
                                duration=mort_vars['Values'].loc['duration'][0]
                                predictor=np.array(mort_vars['Values'].loc['predictor'])
                                response=np.array(mort_vars['Values'].loc['response'])
                                weight=np.ones(len(response))
                                weight[(response<0.1)|(response>0.9)]=10
                                weight[(response<0.02)|(response>0.98)]=500
                                coeffs=self._build_logistic(predictor,response,weight,fit_method='dogbox')
                                mort_dict[names]={'duration':duration,'coeffs':coeffs.tolist()}
                        else:
                            mort={}
                        #Unpack all subdictionaries and combine into master vegparams dictionary for species/community
                        vegparams[i]={**factor, **grow, **duration, **size, **disp, **col, **mort}
                else: 
                    if exten == 'csv':
                        #Add Carra's code here and load into dict called x
                        pass
                    else:
                        raise ValueError('File extension not recognized')
               
                self.vegparams=vegparams
                

                    #Calculate derived size parameters
                    #rsratio=[]
                    #stem_mass_unit=[]
                    #stem_dens_max=[]
                    #ag_mass_max=[]
                    #bg_mass_max=[]
                    #for j in range(len(coms)):
                    #    rsratio.append([self.sizeparams['allocate'][j][0]/(self.sizeparams['allocate'][j][1]+self.sizeparams['allocate'][j][2])])
                    #    stem_mass_unit.append([self.sizeparams['stem_mass_max'][j][0]/self.sizeparams['stem_ht_max'][j][0]])
                    #    stem_dens_max.append([self.sizeparams['pl_dens_max'][j][0]*self.sizeparams['pl_stems_max'][j][0]])
                        #ag_mass_max.append([t1*t2*(1+(t3/t4)) for [t1,t2,[t3,t4,t5]] in 
                        #                    [stem_dens_max[j][0],self.sizeparams['stem_mass_max'][j][0],
                        #                    self.sizeparams['allocate']]])
                        #bg_mass_max.append([t1*t2 for [t1,t2] in [ag_mass_max[j][0],rsratio[j][0]]])
                    #self.sizeparams['rsratio']=rsratio
                    #self.sizeparams['stem_mass_unit']=stem_mass_unit
                    #self.sizeparams['stem_dens_max']=stem_dens_max
                    #self.sizeparams['ag_mass_max']=ag_mass_max
                    #self.sizeparams['bg_mass_max']=bg_mass_max
                
        #Save vegparams dictionary to yaml for future use
        with open(outfile,'w') as outfile:
            yaml.dump(self.vegparams, outfile, default_flow_style=True)

#Private method to read in dataframe from Excel and return a formatted dictionary
    def _makedict(self, df, keys, name):
        #Only read in variables with defined key names
        #temp=df[df['Variable Name'].isin(keys)]
        temp=df.loc[keys]
        if temp['Values'].isnull().values.any():
            msg='Cannot build dictionary for {}. One of the variable values is missing. Please check the input file and try again.'.format(name)
            raise ValueError(msg)
        else:
            #Aggregate and group variables with multiple values as a list
            group_df=temp.groupby(['Variable Name','Descriptor']).agg(pd.Series.tolist)
            #Replace single element lists with float, int, or string
            temp_dict={}
            var_names=group_df.index.levels[0]
            for var in var_names:
                temp_dict[var]={}
                group_df_xs=group_df.xs(var)
                sub_dict_vars=group_df_xs.index.values
                for sub_var in sub_dict_vars:
                    dict_entry=(lambda x: x[0] if len(x)==1 else x)(group_df_xs['Values'].loc[sub_var])
                    temp_dict[var][sub_var]=dict_entry
            return temp_dict

    def replace_part_list_with_dict(self, val_list):
        part_list=['root','leaf','stem','storage','reproductive']
        biomass_list=['root_biomass', 'leaf_biomass','stem_biomass','storage_biomass','repro_biomass']
        replace_dict={}
        for i, part in enumerate(part_list):
            replace_dict[part]=val_list[i]
            #replace_dict[biomass_list[i]]=val_list[i]
        return replace_dict
#Private method to build logistic mortality function for up to five acute mortality factors
    def _build_logistic(self, xs, ys, weights, fit_method):
        S=[]
        ys[ys<=0]=0.0001
        ys[ys>=1]=0.999
        #Check to make sure input arrays are same length
        if len(xs)!=len(ys): 
            msg=('Predictor and response variable arrays must be same length')
            S=[]
            raise ValueError(msg)
        elif len(xs)<=1:
            print(len(xs))
            msg=('Not enough points to generate logistic function. Assuming zero mortality.')
            S=[0,0]
        #Direct solve for coefficients if only two points provided to prevent solver errors
        elif len(xs)==2:
            #Solve for constant (see function below)
            b=-(np.log((1-ys[1])/ys[1])-np.log((1-ys[0])/ys[0]))/(xs[1]-xs[0])
            a=((1-ys[1])/ys[1])/np.exp(-xs[1]*b)
            S=[a,b]
        else:
            #Find first guess a,b by finding points closest to 0.5 survival and near sigmoid limit
            x=[0,0]
            y=[0,0]
            idx_05=(np.abs(ys-0.5)).argmin()
            idx_limit=np.gradient(np.gradient(ys,xs),xs).argmax()
            x[0]=xs[min(idx_05,idx_limit)]
            y[0]=ys[min(idx_05,idx_limit)]
            x[1]=xs[max(idx_05,idx_limit)]
            y[1]=ys[max(idx_05,idx_limit)]
            b_guess=-(np.log((1-ys[1])/ys[1])-np.log((1-ys[0])/ys[0]))/(xs[1]-xs[0])
            a_guess=((1-ys[1])/ys[1])/np.exp(-xs[1]*b_guess)
            S, pcov = curve_fit(self._cfunc, xs, ys, p0=[a_guess, b_guess], sigma=weights, method=fit_method)
        return S
        
    def _cfunc(self, x,a,b):
        return 1/(1+a*np.exp(-b*x))
