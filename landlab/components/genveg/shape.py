
import numpy as np
#Growth form classes and selection method
class Climbing(object):
    def __init__(self):
        pass

class Conical(object):
    def __init__(self):
        pass

    def calculate_crown_volume(self, plants):
        volume=np.pi/12*plants['shoot_sys_width']**2*plants['shoot_height']
        return volume

class Decumbent(object):
    def __init__(self):
        pass

    def calculate_crown_volume(self, plants):
        volume=np.pi/3*plants['shoot_sys_width']**2*plants['shoot_height']
        return volume

class Erect(object):
    def __init__(self):
        pass

    def calculate_crown_volume(self, plants):
        volume=np.pi/4*plants['shoot_sys_width']**2*plants['shoot_height']
        return volume

class Irregular(object):
    def __init__(self):
        pass

class Oval(object):
    def __init__(self):
        pass
    
    def calculate_crown_volume(self, plants):
        volume=np.pi/6*plants['shoot_sys_width']**2*plants['shoot_height']
        return volume

class Prostrate(object):
    def __init__(self):
        pass

class Rounded(object):
    def __init__(self):
        pass

    def calculate_crown_volume(self, plants):
        volume=np.pi/6*plants['shoot_sys_width']**3
        return volume

class Semierect(object):
    def __init__(self):
        pass

class Vase(object):
    def __init__(self):
        pass