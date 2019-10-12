

# fix redirects
import os

import fileinput

_SRC = "source"

with open("build/linkcheck/output.txt", "r") as f:
    link_out = f.readlines()

#%%
for link in link_out:
    if "redirected" in link:
        file_name = os.path.join(_SRC, link.split(":")[0])

        assert os.path.exists(file_name)

        old, new = link.strip().split("] ")[-1].split(" to ")

        with fileinput.FileInput(file_name, inplace=True) as file:
            for line in file:
                line.replace(old, new)
