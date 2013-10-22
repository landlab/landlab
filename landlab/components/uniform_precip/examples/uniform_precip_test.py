from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution


def main():
    print 'We are going to use TrialRun as our class instance.'
    TrialRun = PrecipitationDistribution()
    print 'TrialRun = PrecipitationDistribution()'
    print '\n'

    print "TrialRun's values before initiation..."
    print "Storm Duration is: ", TrialRun.storm_duration, 'hours.'
    print "Interstorm Duration is: ", TrialRun.interstorm_duration, 'hours.'
    print "Storm Depth is: ", TrialRun.storm_depth, 'mm.'
    print "Intensity is: ", TrialRun.intensity, 'mm/hr.' 
    print '\n'
    print 'We should initialize TrialRun... TrialRun.initialize()'
    TrialRun.initialize()
    print '\n'
    print 'What are the mean values read in from the input file?'
    print "Mean Storm Duration is: ", TrialRun.mean_storm, 'hours.'
    print "Interstorm Duration is: ", TrialRun.mean_interstorm, 'hours.'
    print "Storm Depth is: ", TrialRun.mean_storm_depth, 'mm.'
    print "Intensity is: ", TrialRun.mean_intensity, 'mm/hr.'
    print '\n'
    print "Let's see what what the class members are after the first initialization..."
    print "Storm Duration is: ", TrialRun.storm_duration, 'hours.'
    print "Interstorm Duration is: ", TrialRun.interstorm_duration, 'hours.'
    print "Storm Depth is: ", TrialRun.storm_depth, 'mm.'
    print "Intensity is: ", TrialRun.intensity, 'mm/hr.'

    print '\n'
    print 'Now we will update these values using TrialRun.update()'
    TrialRun.update()
    print "Storm Duration is: ", TrialRun.storm_duration, 'hours.'
    print "Interstorm Duration is: ", TrialRun.interstorm_duration, 'hours.'
    print "Storm Depth is: ", TrialRun.storm_depth, 'mm.'
    print "Intensity is: ", TrialRun.intensity, 'mm.'
    
    print '\n'
    print 'Now we are going to generate a time series:'
    TrialRun.get_storm_time_series()
    print TrialRun.storm_time_series


if __name__ == '__main__':
    main()

