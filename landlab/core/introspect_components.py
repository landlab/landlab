import os.path as path
from os import walk
import fnmatch
import glob
import landlab.components as comp
import dircache
import pkgutil
from copy import copy

abspath = path.abspath(comp.__path__[0])
poss_comp_files = []
for root, dirnames, filenames in walk(abspath):
    for filename in fnmatch.filter(filenames, '*.py'):
        poss_comp_files.append(path.join(root, filename))

props_to_strip = [' _name',
                  ' _input_var_names',
                  ' _output_var_names',
                  ' _var_units',
                  ' _var_mapping',
                  ' _var_doc']

comp_elements = {}
problematic_components = set()  # components lacking all the info
last_name = None

for LLcomp in poss_comp_files:
    # print LLcomp
    found_a_name = False
    for prop in props_to_strip:
        lines_captured = []
        start_write = False
        with open(LLcomp, 'r') as inFile:
            for line in inFile:
                if prop in line:
                    start_write = True
                    if prop != ' _name':
                        assert ('{' in line) or ('[' in line)
                    else:
                        found_a_name = True
                if start_write:
                    if not found_a_name:
                        problematic_components.add(LLcomp)
                    nowhite = line.rstrip()
                    nowhite = nowhite.lstrip()
                    no_nl = nowhite.replace('\\', '')
                    lines_captured.append(str(no_nl))
                    if ('}' in line) or (']' in line) or (prop == ' _name'):
                        break
        cat_lines = ''
        # print(lines_captured)
        for expr in lines_captured:
            cat_lines += expr
        cat_lines = cat_lines.replace(" ", "")
        if cat_lines and found_a_name:
            print('EXEC: ', LLcomp)
            exec(cat_lines)  # eval(prop) is now an obj
            if prop is ' _name':
                last_name = eval(prop.lstrip())
                comp_elements[eval(prop.lstrip())] = {}
            else:
                comp_elements[last_name][prop.lstrip()] = copy(eval(prop.lstrip()))

# now, comp_elements is a dict of dict of dicts, where the value at the end
# is the dict or set produced by that component property
