"""
Load veg veg_params
May need to change this to output a landlab configured input file
Want to retain formatted Excel file due to type of data required
"""
import numpy as np
import pandas as pd
import pathlib
import yaml

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
    processes: a list of vegetation processes to initialize parameters for 'plantsize','dispersal','mortality','colonization','storage'. 
               If blank, model will run basic growth only.
    outfile: optional input to allow for custom file name
    """
    def __init__(
        self,
        fpath='None',
        outfile='veg_params.yml',
        processes=[]
    ):
        if fpath=='None':
            self.plant_factors={
                'name':['Corn'],
                'type': [1],
                'ptype': ['C4'],
            }
            self.growthparams= {
                'gs_start': [91],
                'gs_end': [290],
                'gs_sen': [228],
                'res_co': [0.015,0.015,0.03],
                'glu_req': [1.444,1.513,1.463],
                'le_k':[0.02],
                'hi':[9],
                'p_max':[0.055],
                'allocate':[0.2,0.3,0.5],
            }
            if 'plantsize' in processes:
                self.sizeparams={
                    'pl_dens_max': [1],
                    'stems_max': [3],
                    'stem_ht_max': [2.5],
                    'stem_mass_max': [72],
                    'allocate': [0.45,0.3,0.25],
                    'stem_tca': [0.231],
                    'min_mass':[0.01,0.1,0.5]
                }
                if 'dispersal' in processes:
                    self.dispparams={
                        'disp_dist': [2],
                        'disp_size_rat': [0.5],
                        'disp_cost': [0]
                    }
                if 'storage' in processes:
                    self.storparams={
                        'r_wint_die': [0.25],
                        'r_wint_stor': [0.25]
                    }
            if 'colonize' in processes:
                self.colparams={
                    'col_prob': [0.01],
                    'col_dt': [365]
                }
            if 'mortality' in processes:
                self.mortparams={
                    's1_name': 'Mortality factor',
                    's1_days': [365],
                    's1_pred': [1,2,3,4],
                    's1_rate': [0,0.1,0.9,1],
                    's1_weight':[1000,1,1,1000]
                }
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
                        df_in=xlin.parse(i, usecols='B,D',skiprows=[1,4,18,30,34,37,40])
                        keyval=i
                        growth_keys=[
                            'name',
                            'type',
                            'ptype',
                            'gs_start',
                            'gs_end',
                            'senes',
                            'res_co',
                            'glu_req',
                            'le_k',
                            'hi',
                            'p_max',
                            'allocate'
                        ]
                        grow=df_in[df_in['Variable Name'].isin(growth_keys)]
                        group_grow=grow.groupby('Variable Name').agg(pd.Series.tolist)
                        group_grow.rename(columns={'Values': 'growparams'}, inplace=True)
                        grow=group_grow.to_dict()
                        if 'plantsize' in processes:
                            size_keys=[
                                'pl_dens_max',
                                'pl_stems_max',
                                'stem_ht_max',
                                'stem_mass_max',
                                'stem_tca',
                                'min_mass'
                            ]
                            size=df_in[df_in["Variable Name"].isin(size_keys)]
                            group_size=size.groupby('Variable Name').agg(pd.Series.tolist)
                            group_size.rename(columns={'Values': 'sizeparams'}, inplace=True)
                            size=group_size.to_dict()
                            if 'dispersal' in processes:
                                disp_keys=[
                                    'disp_dist',
                                    'disp_size_rat',
                                    'disp_cost'
                                ]
                                disp=df_in[df_in["Variable Name"].isin(disp_keys)]
                                group_disp=disp.groupby('Variable Name').agg(pd.Series.tolist)
                                group_disp.rename(columns={'Values': 'dispparams'}, inplace=True)
                                disp=group_disp.to_dict()
                            else:
                                disp={{}}
                            if 'storage' in processes:
                                stor_keys=[
                                    'r_wint_die',
                                    'r_wint_stor'
                                ]
                                stor=df_in[df_in["Variable Name"].isin(stor_keys)]
                                group_stor=stor.groupby('Variable Name').agg(pd.Series.tolist)
                                group_stor.rename(columns={'Values': 'storparams'}, inplace=True)
                                stor=group_stor.to_dict()
                            else:
                                stor={{}}
                        else:
                            size={{}}
                        if 'colonize' in processes:
                            col_keys=[
                                'col_prob',
                                'col_dt'
                            ]
                            col=df_in[df_in["Variable Name"].isin(col_keys)]
                            group_col=col.groupby('Variable Name').agg(pd.Series.tolist)
                            group_col.rename(columns={'Values': 'colparams'}, inplace=True)
                            col=group_col.to_dict()
                        else:
                            col={{}}
                        if 'mortality' in processes:
                            mort_keys=[
                                's1_name',
                                's1_days',
                                's1_pred',
                                's1_rate',
                                's1_weight'
                                's2_name',
                                's2_days',
                                's2_pred',
                                's2_rate',
                                's2_weight',
                                's3_name',
                                's3_days',
                                's3_pred',
                                's3_rate',
                                's3_weight',
                                's4_name',
                                's4_days',
                                's4_pred',
                                's4_rate',
                                's4_weight',
                                's5_name',
                                's5_days',
                                's5_pred',
                                's5_rate',
                                's5_weight'
                            ]
                            mort=df_in[df_in["Variable Name"].isin(mort_keys)]
                            group_mort=mort.groupby('Variable Name').agg(pd.Series.tolist)
                            group_mort.rename(columns={'Values': 'mortparams'}, inplace=True)
                            mort=group_mort.to_dict()
                        else:
                            mort={{}}
                        vegparams={keyval:{**grow,**size,**disp,**stor,**col,**mort}}
                else: 
                    if exten == 'csv':
                        #Add Carra's code here and load into dict called x
                        pass
                    else:
                        raise ValueError('File extension not recognized')
                
                #Make list of vegparams keys to create new dictionary    
               
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
                

        with open(outfile,'w') as outfile:
            yaml.dump(self.vegparams, outfile, default_flow_style=True)            