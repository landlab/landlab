import os

from landlab import load_params

required_params = ['number_of_node_columns','dx','run_duration','dt','output_interval']

DEFAULT_PARAMS = load_params('default.yaml')


def Yaml_File_Check(file):
    if 'yaml' in file:
        from landlab import load_params
        planform = load_params(file)
        for element in required_params:
            missing_required_params = []
            if element not in list(planform.keys()):
                 missing_required_params.append(element)
            for value in missing_required_params:
                planform[value] = DEFAULT_PARAMS[value] 
    else:
        planform = DEFAULT_PARAMS

    return planform


# An alternative that does the same thing.
def yaml_file_check(file):
    root, ext = os.path.splitext(file)
    if ext == '.yaml':
        planform = load_params(file)
        for element in required_params:
            planform.setdefault(element, DEFAULT_PARAMS[element])
    else:
        planform = DEFAULT_PARAMS

    return planform