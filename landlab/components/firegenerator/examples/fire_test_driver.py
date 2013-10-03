from landlab.components.firegenerator.fire_generator import FireGenerator
from matplotlib import pyplot as plt


FG = FireGenerator()
FG.initialize()
FG.get_scale_parameter()
FG.generate_fire_time_series()
FG.fires.sort()
yaxis = range(FG.total_run_time)
plt.plot(FG.fires, yaxis)
plt.figure(2)
plt.hist(FG.fires)
plt.show()
