from .dispersal import *

#Growth form classes and selection method
class Bunch(Seed):
    def __init__(self):
        pass
    
    def branch(self):
        print('Limited lateral branching due to clumping')
    
    def disperse(self, plants):
        return plants

class Colonizing(Random):
    def __init__(self):
        pass

    def branch(self):
        print('No branching annual')
    
    def disperse(self, plants):
        return plants

class Multiplestems(Seed):
    def __init__(self):
        pass

    def branch(self):
        print('Create two or more main stems at or near soil surface')
    
    def disperse(self, plants):
        return plants

class Rhizomatous(Clonal):
    def __init__(self):
        pass
    def branch(self):
        print('Tiller via rhizomes')

    def disperse(self, plants):
        return plants

class Singlecrown(Seed):
    def __init__(self):
        pass

    def branch(self):
        print('Herbaceous plant with one persistent base')

    def disperse(self, plants):
        return plants

class Singlestem(Seed):
    def __init__(self):
        pass

    def branch(self):
        print('Plant develops one stem like a tree or a corn plant')

    def disperse(self, plants):
        return plants

class Stoloniferous(Clonal):
    def __init__(self):
        pass

    def branch(self):
        print('Branches via stolons')

    def disperse(self, plants):
        return plants

class Thicketforming(Seed):
    def __init__(self):
        pass

    def branch(self):
        print('Limited lateral branching due to dense thickets')

    def disperse(self, plants):
        return plants