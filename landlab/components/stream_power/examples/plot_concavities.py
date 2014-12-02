from pylab import figure, plot, xlabel, ylabel, title

limit_of_steepening=[]
for i in xrange(len(A_profiles)):
    counter=0
    steepest=0.
    for n in S_profiles[i]:
        #print n, steepest
        counter+=1
        if n<steepest-0.001*np.amax(S_profiles):
            limit_of_steepening.append(counter)
            break
        steepest=n

print limit_of_steepening

concavities = []
convexity_x = []
max_x = []
for i in xrange(len(limit_of_steepening)):
    lim=limit_of_steepening[i]
    polyvals = np.polyfit(np.log10(A_profiles[i][1:lim]), np.log10(S_profiles[i][1:lim]), 1)
    if polyvals[0]<0:
        concavities.append(-polyvals[0])
        convexity_x.append(x_profiles[i][0][lim-1])
        max_x.append(np.amax(x_profiles[0]))

concavities = np.array(concavities)
convexity_x = np.array(convexity_x)
max_x = np.array(max_x)

figure('long_profiles')
for i in xrange(0,67):
    loglog(A_profiles[i],S_profiles[i], 'x-')

yunnan_propx = np.loadtxt('yunnan_proplength.txt')
yunnan_theta = np.loadtxt('yunnan_theta.txt')
fagaras_propx = np.loadtxt('fagaras_proplength.txt')
fagaras_theta = np.loadtxt('fagaras_theta.txt')
ladakh_propx = np.loadtxt('ladakh_proplength.txt')
ladakh_theta = np.loadtxt('ladakh_theta.txt')
#this data is u=0.0001->0.0005, no threshold, explicit, small random K var from 1e-6
pureDL_nothresh_propx = np.loadtxt('pureDLnothresh_proplength.txt')
pureDL_nothresh_theta = np.loadtxt('pureDLnothresh_theta.txt')
seddepNMG_propx = np.loadtxt('seddepNMG_proplength.txt')
seddepNMG_theta = np.loadtxt('seddepNMG_theta.txt')

figure('concavities')
plot(yunnan_propx, np.log(yunnan_theta), 'o')
plot(fagaras_propx, np.log(fagaras_theta), 'o')
plot(ladakh_propx, np.log(ladakh_theta), 'o')
plot(1.-convexity_x/max_x, np.log(concavities), '*-')
plot(seddepNMG_propx, np.log(seddepNMG_theta), 'x-')
plot(pureDL_nothresh_propx, np.log(pureDL_nothresh_theta), '.-')
xlabel('relative position of knickzone down channel')
ylabel('natural log of concavity downstream of knickzone')
title('Elevated concavities due to propagating knickzone in sed flux dependent incision')