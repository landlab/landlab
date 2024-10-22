import numpy as np
import pytest


@pytest.fixture(scope="session", autouse=True)
def set_numpy_printoptions():
    np.set_printoptions(legacy="1.25")
