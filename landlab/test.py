from decorator import implements

class Bmi(object):
    input_var_names = []

    def initialize(self):
        pass
    def run(self):
        pass

@implements(Bmi)
class MyBmi(object):
    input_var_names = [
        'var_string',
    ]

    def initialize(self):
        return 'Done'
    def run(self):
        return until

print MyBmi.__implements__
