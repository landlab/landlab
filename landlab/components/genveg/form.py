from .dispersal import *

#Growth form classes and selection method
class Bunch(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Limited lateral branching due to clumping')

class Colonizing(Random):
    def __init__(self):
        pass
    def branch(self):
        print('No branching annual')

class Multiplestems(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Create two or more main stems at or near soil surface')

class Rhizomatous(Clonal):
    def __init__(self):
        pass
    def branch(self):
        print('Tiller via rhizomes')

class Singlecrown(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Herbaceous plant with one persistent base')

class Singlestem(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Plant develops one stem like a tree or a corn plant')

class Stoloniferous(Clonal):
    def __init__(self):
        pass
    def branch(self):
        print('Branches via stolons')

class Thicketforming(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Limited lateral branching due to dense thickets')