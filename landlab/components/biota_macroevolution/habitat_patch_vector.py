"""BiotaEvolver HabitatPatchVector object.
"""


class HabitatPatchVector(object):

    # Define cardinality cases.
    NONE_TO_ONE = 'none-to-one'
    NONE_TO_MANY = 'none-to-many'
    ONE_TO_NONE = 'one-to-none'
    ONE_TO_ONE = 'one-to-one'
    ONE_TO_MANY = 'one-to-many'
    MANY_TO_NONE = 'many-to-none'
    MANY_TO_ONE = 'many-to-one'
    MANY_TO_MANY = 'many-to-many'

    def __init__(self, prior_patch_count, new_patch_count):

        self.origin = None
        self.destination = []

        # Set cardinality.

        if prior_patch_count == 0:
            if new_patch_count == 1:
                self.cardinality = self.NONE_TO_ONE
            elif new_patch_count > 1:
                self.cardinality = self.NONE_TO_MANY

        if prior_patch_count == 1:
            if new_patch_count == 0:
                self.cardinality = self.ONE_TO_NONE
            elif new_patch_count == 1:
                self.cardinality = self.ONE_TO_ONE
            elif new_patch_count > 1:
                self.cardinality = self.ONE_TO_MANY

        elif prior_patch_count > 1:
            if new_patch_count == 0:
                self.cardinality = self.MANY_TO_NONE
            elif new_patch_count == 1:
                self.cardinality = self.MANY_TO_ONE
            elif new_patch_count > 1:
                self.cardinality = self.MANY_TO_MANY
