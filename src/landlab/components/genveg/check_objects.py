import numpy as np


class UnitTestChecks:
    def is_negative_present(self, test_value, print_variable):
        if isinstance(test_value, np.ndarray):
            if test_value[test_value < 0].size > 0:
                raise ValueError(
                    f"A negative value was found in {print_variable} array"
                )
        else:
            if test_value < 0:
                raise ValueError(f"{print_variable} is negative")

    def is_zero(self, test_value, print_variable):
        if test_value == 0:
            raise ValueError(f"{print_variable} is zero, cannot divide by zero")
