import os

from hypothesis import settings

settings.register_profile("ci", deadline=None)
settings.load_profile(os.getenv("HYPOTHESIS_PROFILE", "default"))
