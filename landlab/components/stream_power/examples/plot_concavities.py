from pylab import figure, plot, xlabel, ylabel, title, loglog, show
import numpy as np

limit_of_steepening=[]
max_S = 0.
for i in xrange(len(S_profiles)):
    if np.amax(S_profiles[i])>max_S:
        max_S = np.amax(S_profiles[i])

for i in xrange(len(S_profiles)):
    counter=0
    steepest=0.
    for n in S_profiles[i]:
        #print n, steepest
        counter+=1
        if n<steepest-0.002*max_S:#np.amax(S_profiles[i]):
        #if n<0.7*steepest:
            limit_of_steepening.append(counter)
            break
        steepest=n
start_node = 1

#if using a manually picked dictionary:
#copy this chunk out if not
#limit_of_steepening = [0]*len(A_profiles)
#for i in kz_pts_dict.keys():
#    limit_of_steepening[i] = kz_pts_dict[i]
#start_node = 5 #because there are awful edge effects if using storms

print limit_of_steepening

concavities = []
convexity_x = []
max_x = []
for i in xrange(len(limit_of_steepening)):
    lim=limit_of_steepening[i]
    if lim != 0:
        polyvals = np.polyfit(np.log10(A_profiles[i][start_node:lim]), np.log10(S_profiles[i][start_node:lim]), 1)
        if polyvals[0]<0:
            concavities.append(-polyvals[0])
            convexity_x.append(x_profiles[i][lim-1])
            max_x.append(np.amax(x_profiles[i]))

concavities = np.array(concavities)
convexity_x = np.array(convexity_x)
max_x = np.array(max_x)

figure('long_profiles')
for i in xrange(len(A_profiles)):
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
stormsDL_thresh_propx = np.loadtxt('stormsDLthresh1e-1_proplength.txt')
stormsDL_thresh_theta = np.loadtxt('stormsDLthresh1e-1_theta.txt')
stormsDL_nothresh_propx = np.loadtxt('stormsDLnothresh_proplength.txt')
stormsDL_nothresh_theta = np.loadtxt('stormsDLnothresh_theta.txt')
bevel1seddep_propx = np.loadtxt('seddepNMG0.0001bevel_proplength.txt')
bevel1seddep_theta = np.loadtxt('seddepNMG0.0001bevel_theta.txt')
bevel2seddep_propx = np.loadtxt('seddepNMG0.0006bevel_proplength.txt') #no perturbation... we're not crossing the hump
bevel2seddep_theta = np.loadtxt('seddepNMG0.0006bevel_theta.txt')
aparabolicNMG2x10_propx = np.loadtxt('aparabolicNMG20001x10_proplength.txt')
aparabolicNMG2x10_theta = np.loadtxt('aparabolicNMG20001x10_theta.txt')
aparabolicNMG2x5_propx = np.loadtxt('aparabolicNMG20001x5_proplength.txt')
aparabolicNMG2x5_theta = np.loadtxt('aparabolicNMG20001x5_theta.txt')
aparabolicNMG4x100_propx = np.loadtxt('aparabolicNMG4000001x100_proplength.txt')
aparabolicNMG4x100_theta = np.loadtxt('aparabolicNMG4000001x100_theta.txt')
aparabolicNMG4x10_propx = np.loadtxt('aparabolicNMG4000001x10_proplength.txt')
aparabolicNMG4x10_theta = np.loadtxt('aparabolicNMG4000001x10_theta.txt')
aparabolicNMG5x10_propx = np.loadtxt('aparabolicNMG500001x10_proplength.txt')
aparabolicNMG5x10_theta = np.loadtxt('aparabolicNMG500001x10_theta.txt')

figure('concavities')
plot(yunnan_propx, np.log(yunnan_theta), 'o')
plot(fagaras_propx, np.log(fagaras_theta), 'o')
plot(ladakh_propx, np.log(ladakh_theta), 'o')
plot(seddepNMG_propx, np.log(seddepNMG_theta), 'x-', label='U=0.0001,5*dU')
plot(bevel1seddep_propx, np.log(bevel1seddep_theta), 'x-', label='bevel_lowangle')
plot(aparabolicNMG2x10_propx[1:], np.log(aparabolicNMG2x10_theta[1:]), 'x-', label='10*U,2*K,10*dU')
plot(aparabolicNMG2x5_propx, np.log(aparabolicNMG2x5_theta), 'x-', label='10*U,2*K,5*dU')
plot(aparabolicNMG4x100_propx, np.log(aparabolicNMG4x100_theta), 'x-', label='0.1*U,100*dU')
plot(aparabolicNMG4x10_propx, np.log(aparabolicNMG4x10_theta), 'x-', label='0.1*U,10*dU')
plot(aparabolicNMG5x10_propx, np.log(aparabolicNMG5x10_theta), 'x-', label='0.2*K,10*dU')
plot(bevel2seddep_propx[1:], np.log(bevel2seddep_theta[1:]), 'x-', label='bevel_highangle') #first val is a bad pick
plot(pureDL_nothresh_propx, np.log(pureDL_nothresh_theta), '.-', label='pureDL')
plot(stormsDL_thresh_propx, np.log(stormsDL_thresh_theta), '+-', label='pureDL_storms_threshold')
plot(stormsDL_nothresh_propx, np.log(stormsDL_nothresh_theta), '+-', label='pureDL_storms')
plot(1.-convexity_x/max_x, np.log(concavities), '*-', label='active_run')
xlabel('relative position of knickzone down channel')
ylabel('natural log of concavity downstream of knickzone')
title('Concavities donstream of propagating knickzones')
pylab.legend(loc=2)

figure('concavities_nolog')
plot(yunnan_propx, yunnan_theta, 'o')
plot(fagaras_propx, fagaras_theta, 'o')
plot(ladakh_propx, ladakh_theta, 'o')
plot(seddepNMG_propx, seddepNMG_theta, 'x-')
plot(bevel1seddep_propx, bevel1seddep_theta, 'x-')
plot(aparabolicNMG2x10_propx[1:], aparabolicNMG2x10_theta[1:], 'x-')
plot(aparabolicNMG2x5_propx, aparabolicNMG2x5_theta, 'x-')
plot(aparabolicNMG4x100_propx, aparabolicNMG4x100_theta, 'x-')
plot(aparabolicNMG4x10_propx, aparabolicNMG4x10_theta, 'x-')
plot(aparabolicNMG5x10_propx, aparabolicNMG5x10_theta, 'x-')
plot(bevel2seddep_propx[1:], bevel2seddep_theta[1:], 'x-') #first val is a bad pick
plot(pureDL_nothresh_propx, pureDL_nothresh_theta, '.-')
plot(stormsDL_thresh_propx, stormsDL_thresh_theta, '+-')
plot(stormsDL_nothresh_propx, stormsDL_nothresh_theta, '+-')
plot(1.-convexity_x/max_x, concavities, '*-')
xlabel('relative position of knickzone down channel')
ylabel('concavity downstream of knickzone')
title('Concavities donstream of propagating knickzones')