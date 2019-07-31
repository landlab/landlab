import numpy as np
from pylab import figure, gca, legend, plot, xlim, ylim

yunnan_propx = np.loadtxt("yunnan_proplength.txt")
yunnan_theta = np.loadtxt("yunnan_theta.txt")
fagaras_propx = np.loadtxt("fagaras_proplength.txt")
fagaras_theta = np.loadtxt("fagaras_theta.txt")
ladakh_propx = np.loadtxt("ladakh_proplength.txt")
ladakh_theta = np.loadtxt("ladakh_theta.txt")
# this data is u=0.0001->0.0005, no threshold, explicit, small random K var from 1e-6
pureDL_nothresh_propx = np.loadtxt("pureDLnothresh_proplength.txt")
pureDL_nothresh_theta = np.loadtxt("pureDLnothresh_theta.txt")
seddepNMG_propx = np.loadtxt("seddepNMG_proplength.txt")
seddepNMG_theta = np.loadtxt("seddepNMG_theta.txt")
stormsDL_thresh_propx = np.loadtxt("stormsDLthresh1e-1_proplength.txt")
stormsDL_thresh_theta = np.loadtxt("stormsDLthresh1e-1_theta.txt")
stormsDL_nothresh_propx = np.loadtxt("stormsDLnothresh_proplength.txt")
stormsDL_nothresh_theta = np.loadtxt("stormsDLnothresh_theta.txt")
bevel1seddep_propx = np.loadtxt("seddepNMG0.0001bevel_proplength.txt")
bevel1seddep_theta = np.loadtxt("seddepNMG0.0001bevel_theta.txt")
bevel2seddep_propx = np.loadtxt(
    "seddepNMG0.0006bevel_proplength.txt"
)  # no perturbation... we're not crossing the hump
bevel2seddep_theta = np.loadtxt("seddepNMG0.0006bevel_theta.txt")
aparabolicNMG2x10_propx = np.loadtxt("aparabolicNMG20001x10_proplength.txt")
aparabolicNMG2x10_theta = np.loadtxt("aparabolicNMG20001x10_theta.txt")
aparabolicNMG2x5_propx = np.loadtxt("aparabolicNMG20001x5_proplength.txt")
aparabolicNMG2x5_theta = np.loadtxt("aparabolicNMG20001x5_theta.txt")
aparabolicNMG4x100_propx = np.loadtxt("aparabolicNMG4000001x100_proplength.txt")
aparabolicNMG4x100_theta = np.loadtxt("aparabolicNMG4000001x100_theta.txt")
aparabolicNMG4x10_propx = np.loadtxt("aparabolicNMG4000001x10_proplength.txt")
aparabolicNMG4x10_theta = np.loadtxt("aparabolicNMG4000001x10_theta.txt")
aparabolicNMG5x10_propx = np.loadtxt("aparabolicNMG500001x10_proplength.txt")
aparabolicNMG5x10_theta = np.loadtxt("aparabolicNMG500001x10_theta.txt")

figure("concavities_just_one_sde")
plot(seddepNMG_propx, seddepNMG_theta, "rx-", label="sediment flux dependent")
gca().set_yscale("log")
xlim([0, 1])
y_scale = gca().get_ylim()
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")

figure("concavities_just_one_DL")
plot(
    pureDL_nothresh_propx, pureDL_nothresh_theta, "b.-", label="pure detachment limited"
)
gca().set_yscale("log")
xlim([0, 1])
ylim(y_scale)
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")
legend(loc=2)

figure("concavities_one_sde_one_DL")
plot(
    pureDL_nothresh_propx, pureDL_nothresh_theta, "b.-", label="pure detachment limited"
)
plot(seddepNMG_propx, seddepNMG_theta, "rx-", label="sediment flux dependent")
gca().set_yscale("log")
xlim([0, 1])
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")
legend(loc=2)

figure("all_models")
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")
plot(pureDL_nothresh_propx, pureDL_nothresh_theta, "b.-", label="pureDL")
plot(stormsDL_nothresh_propx, stormsDL_nothresh_theta, "b+-", label="pureDL_storms")
plot(stormsDL_thresh_propx, stormsDL_thresh_theta, "b*-", label="pureDL_storms_thresh")
plot(seddepNMG_propx, seddepNMG_theta, "rx-", label="U=0.0001,K=5*10^-5,5*dU")
plot(
    aparabolicNMG2x10_propx[1:],
    aparabolicNMG2x10_theta[1:],
    "r+-",
    label="10*U,2*K,10*dU",
)
plot(aparabolicNMG2x5_propx, aparabolicNMG2x5_theta, "rv-", label="10*U,2*K,5*dU")
plot(aparabolicNMG4x100_propx, aparabolicNMG4x100_theta, "r<-", label="0.1*U,100*dU")
plot(aparabolicNMG4x10_propx, aparabolicNMG4x10_theta, "r>-", label="0.1*U,10*dU")
plot(aparabolicNMG5x10_propx, aparabolicNMG5x10_theta, "r^-", label="0.2*K,10*dU")
plot(bevel1seddep_propx, bevel1seddep_theta, "rp-", label="bevel_lowangle")
plot(
    bevel2seddep_propx[1:], bevel2seddep_theta[1:], "rh-", label="bevel_highangle"
)  # first val is a bad pick
gca().set_yscale("log")
y_scale_all = gca().get_ylim()
xlim([0, 1])
legend(loc=2)

figure("all_DL")
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")
plot(pureDL_nothresh_propx, pureDL_nothresh_theta, "b.-", label="pureDL")
plot(stormsDL_nothresh_propx, stormsDL_nothresh_theta, "b+-", label="pureDL_storms")
plot(stormsDL_thresh_propx, stormsDL_thresh_theta, "b*-", label="pureDL_storms_thresh")
gca().set_yscale("log")
xlim([0, 1])
ylim(y_scale_all)
legend(loc=2)

figure("all_model_just_Ladakh")
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")
plot(
    pureDL_nothresh_propx,
    pureDL_nothresh_theta,
    ":",
    color="0.6",
    label="pure detachment limited",
)
plot(stormsDL_nothresh_propx, stormsDL_nothresh_theta, ":", color="0.6")
plot(stormsDL_thresh_propx, stormsDL_thresh_theta, ":", color="0.6")
plot(seddepNMG_propx, seddepNMG_theta, "-", color="0.3", label="sed flux dependent")
plot(aparabolicNMG2x10_propx[1:], aparabolicNMG2x10_theta[1:], "-", color="0.3")
plot(aparabolicNMG2x5_propx, aparabolicNMG2x5_theta, "-", color="0.3")
plot(aparabolicNMG4x100_propx, aparabolicNMG4x100_theta, "-", color="0.3")
plot(aparabolicNMG4x10_propx, aparabolicNMG4x10_theta, "-", color="0.3")
plot(aparabolicNMG5x10_propx, aparabolicNMG5x10_theta, "-", color="0.3")
plot(bevel1seddep_propx, bevel1seddep_theta, "-", color="0.3")
plot(
    bevel2seddep_propx[1:], bevel2seddep_theta[1:], "-", color="0.3"
)  # first val is a bad pick
plot(ladakh_propx, ladakh_theta, "bo", label="Ladakh field data")
gca().set_yscale("log")
xlim([0, 1])
ylim(y_scale_all)
legend(loc=2)

figure("all_model_all_data")
plot(np.array([0., 1.]), 0.5 * np.ones(2.), "k--")
plot(
    pureDL_nothresh_propx,
    pureDL_nothresh_theta,
    ":",
    color="0.6",
    label="pure detachment limited",
)
plot(stormsDL_nothresh_propx, stormsDL_nothresh_theta, ":", color="0.6")
plot(stormsDL_thresh_propx, stormsDL_thresh_theta, ":", color="0.6")
plot(seddepNMG_propx, seddepNMG_theta, "-", color="0.3", label="sed flux dependent")
plot(aparabolicNMG2x10_propx[1:], aparabolicNMG2x10_theta[1:], "-", color="0.3")
plot(aparabolicNMG2x5_propx, aparabolicNMG2x5_theta, "-", color="0.3")
plot(aparabolicNMG4x100_propx, aparabolicNMG4x100_theta, "-", color="0.3")
plot(aparabolicNMG4x10_propx, aparabolicNMG4x10_theta, "-", color="0.3")
plot(aparabolicNMG5x10_propx, aparabolicNMG5x10_theta, "-", color="0.3")
plot(bevel1seddep_propx, bevel1seddep_theta, "-", color="0.3")
plot(
    bevel2seddep_propx[1:], bevel2seddep_theta[1:], "-", color="0.3"
)  # first val is a bad pick
plot(yunnan_propx, yunnan_theta, "sr", label="Red River field data")
plot(ladakh_propx, ladakh_theta, "bo", label="Ladakh field data")
plot(fagaras_propx, fagaras_theta, "vg", label="Fagaras field data")
gca().set_yscale("log")
xlim([0, 1])
ylim(y_scale_all)
legend(loc=2)
