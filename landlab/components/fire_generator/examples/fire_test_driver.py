from landlab.components.fire_generator.generate_fire import FireGenerator
from matplotlib import pyplot as plt

print 'Initializing the FireGenerator() class from fire_generator.py'
FG = FireGenerator()
FG.initialize()

print 'Finding unknown scale parameter given the mean fire recurrence value'
FG.get_scale_parameter()

print 'Finding a time series based on the Weibull distribution'
FG.generate_fire_time_series()

print 'Sorting the fire series to test the distribution...'
FG.fires.sort()

print 'Plotting the histogram to show the distribution.'
yaxis = range(FG.total_run_time)
plt.plot(FG.fires, yaxis)
plt.hist(FG.fires)
plt.show()
