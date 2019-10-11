#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Thu Oct 10 13:03:01 2019.

@author: barnhark
"""

import textwrap
from tabulate import tabulate
import numpy as np
import pandas as pd
from landlab.components import COMPONENTS

#%%
    
out = []
for comp in COMPONENTS:
    for name in comp._info:
        temp = {'component': comp.__name__, 
                'field': name}
        for key in comp._info[name].keys(): 
            temp[key] = comp._info[name][key]
        out.append(temp)
df = pd.DataFrame(out)

unique_fields = df.field.unique().astype(str)
#%%

# Table 1, List field name and description
table = []
headers = ("Field Name", "Definition")
for field in np.sort(unique_fields):
    where = df.field == field
    doc = df.doc[where].values[0]
    
    doc_split = textwrap.wrap(doc, 40)
    table.append((field, "\n".join(doc_split)))

out = tabulate(table, headers, tablefmt="grid")
        


#%%
# Table 2, List field IO by component

table = []
headers = ("Field Name", "Provided By", "Used By")
for field in np.sort(unique_fields):

    where = df.field == field
    sel = df[where]
    in_comps = []
    out_comps = []
    
    for i in range(sel.shape[0]):
        name = sel.component.values[i]
        io = sel.intent.values[i]
        mapping = sel.mapping.values[i]
        
        listing = ":py:class:`~landlab.components."+name + "` ("+mapping+")"
        
        if "in" in io:
            in_comps.append(listing)
        if "out" in io:
            out_comps.append(listing)
      
            
    doc_split = textwrap.wrap(doc, 40)
    table.append((field, "\n".join(in_comps), "\n".join(out_comps)))

out = tabulate(table, headers, tablefmt="grid")

with open("temp.txt", "w") as f:
    f.write(out)