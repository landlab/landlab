from landlab.data_record import DataRecord
from xarray import load_mfdataset

# See if eric thinks this is the right place to put this or if it shoudl go in
# landlab.io.

"""
Goal here is to be able to run code that looks something like the following:

# create network model grid, parcels, and NetworkSedimentTransporter

# run model:

Examples
--------
>>> from landlab import NetworkModelGrid
>>> from landlab.components import NetworkSedimentTransporter
>>> from landlab.data_record import DataRecord
>>> from landlab.data_record.io import dump, load_from_multiple_files

TODO: create grid, model, parcels.

# While running be able to dump portion of the parcels to netcdf.

>>> files = []
>>> for t in time:
>>> 	nst.run_one_step(dt)
>>> 	if t%10 == 0: # every ten timesteps write data out.
>>>
>>>       	# write out all but present time to a file, remove prior times from
>>>       	# parcels.
>>>         # define a filename
>>>       	dump(nst.split(), filename)
>>>       	files.append(filename)

Afterwards be able to combine them (including correct join method,
concatenating on the time dimension, etc)

>>> combine_nst_data = load_from_multiple_files(files)

"""


def dump(datarecord, filename):
    """
    """
    pass

def load_from_multiple_files(files):
    """
    """
    # use load mf dataset from xarray.
    pass
