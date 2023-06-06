# Leaf retention classes and selection method
# This may become an attribute of the Duration class unless we define some specific method here
class Evergreen(object):
    def __init__(self):
        self.keep_green_parts = True


class Deciduous(object):
    def __init__(self):
        self.keep_green_parts = False
