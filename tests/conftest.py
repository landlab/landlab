import os

from hypothesis import settings

settings.register_profile("ci", deadline=None, max_examples=10)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))
