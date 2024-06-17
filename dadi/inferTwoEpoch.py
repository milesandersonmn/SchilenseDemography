import nlopt
import dadi
import pickle
import matplotlib.pyplot as plt
import os
import dadi.cuda
dadi.cuda_enabled(False)

os.chdir("/home/milesanderson/PhD/dadi")

def two_epoch_mig(params,ns,pts):
    """
    nu1: size of population 1 after split
    nu2: size of population 2 after split
    Tsplit: time of population split
    mSplit: migration rate after split
    nu1B: size of population 1 during bottleneck
    nu2B: size of population 2 during bottleneck
    TB: time of bottleneck in population
    mB: migration rate during bottleneck
    """
    nu1, nu2, Tsplit, mSplit, nu1B, nu2B, TB, mB=params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Tsplit, nu1=nu1, nu2=nu2, m12=mSplit, m21=mSplit)  # sizes after split
    nu1_func = lambda t: nu1 * (nu1B/nu1) ** (t/TB)
    nu2_func = lambda t: nu2 * (nu2B/nu2) ** (t/TB)
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1_func, nu2=nu2_func, m12= mB, m21= mB)  # simultaneous bottleneck
 
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


# If we have multiple datasets to work with, we can make a variable to store the name of the dataset that we can easily change and redo the inference with a different dataset.
dataset = 'CV2xCV3'

# Load the synonymous frequency spectrum
data_fs = dadi.Spectrum.from_file('/home/milesanderson/PhD/dadi/data/fs/'+dataset+'.fs')

# Retrive the sample sizes from the data
ns = data_fs.sample_sizes

# Define the grid points based on the sample size.
# For smaller data (largest sample size is about <100) [ns+20, ns+30, ns+40] is a good starting point.
# For larger data (largest sample size is about >=100) or for heavily down projected data [ns+100, ns+110, ns+120] is a good starting point.
pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

demo_model = two_epoch_mig

# If the data is unfolded (the ancestral allele was known), as the example data is,
# wrap the demographic model in a function that adds a parameter to estimate the rate of misidentification.
demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)

# Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

# Define starting parameters
params = [1, 1, 0.015, 1, 1, 1, 0.015, 1, 0.01]

# Define boundaries of optimization.
# It is a good idea to have boundaries to avoid optimization
# from trying parameter sets that are time consuming without
# nessicarily being correct.
# If optimization infers parameters very close to the boundaries, we should increase them.
lower_bounds = [0.4, 0.05, 1e-3, 1e-3, 0.05, 0.02, 1e-3, 1e-3, 1e-3]
upper_bounds = [5, 5, 0.5, 10, 5, 5, 0.5, 10, 1]

# Create or append to an file to store optimization results
for i in range(100):
    try:
      fid = open('results/'+dataset+'_demo_fits_two_epoch_mig.txt','a')
    except:
      fid = open('results/'+dataset+'_demo_fits_two_epoch_mig.txt','w')


    # Perturb parameters
    # Optimizers dadi uses are mostly deterministic
    # so we will want to randomize parameters for each optimization.
    # It is recommended to optimize at least 20 time if we only have
    # our local machin, 100 times if we have access to an HPC.
    # If we want a single script to do multiple runs, we will want to
    # start a for loop here
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                                  lower_bound=lower_bounds)

    # Run optimization
    # At the end of the optimization we will get the
    # optimal parameters and log-likelihood.
    # We can modify verbose to watch how the optimizer behaves,
    # what number we pass it how many evaluations are done
    # before the evaluation is printed.

    popt, ll_model = dadi.Inference.opt(p0, data_fs, demo_model_ex, pts_l,
                                        lower_bound=lower_bounds,
                                        upper_bound=upper_bounds,
                                        algorithm=nlopt.LN_BOBYQA,
                                        maxeval=10000, verbose=100)

    # Calculate the synonymous theta
    model_fs = demo_model_ex(popt, ns, pts_l)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)

    # Write results to fid
    res = [ll_model] + list(popt) + [theta0]
    fid.write('\t'.join([str(ele) for ele in res])+'\n')
    fid.close()
    
inferredParams = open('results/'+dataset+'_demo_fits_two_epoch_mig.txt')
lines = inferredParams.readlines()
lines.sort()
bestParams = lines[0]
bestParams = bestParams.split('\t')
bestParams = list(map(float, bestParams))
popt = bestParams[1:10]
model_fs = demo_model_ex(popt, ns, pts_l)

import matplotlib.pyplot as plt
fig = plt.figure(219033)
fig.clear()
dadi.Plotting.plot_2d_comp_multinom(model_fs, data_fs)
fig.savefig('results/'+dataset+'_demo_plot_two_epoch_mig.png')

# Load bootstraped frequency spectrum (if they haven't been made, there is an example in the "Creating a frequency spectrum" section)
import glob
boots_fids = glob.glob('data/fs/bootstraps/'+dataset+'/'+dataset+'.boot_*.fs')
boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

# Godambe uncertainties
# Will contain uncertainties for the
# estimated demographic parameters and theta.
# Start a file to contain the confidence intervals
fi = open('results/'+dataset+'_demographic_confidence_intervals_two_epoch.txt','w')
fi.write('Optimized parameters: {0}\n\n'.format(popt))

# we want to try a few different step sizes (eps) to see if
# uncertainties very wildly with changes to step size.
for eps in [0.01, 0.001, 0.0001]:
    uncerts_adj = dadi.Godambe.GIM_uncert(demo_model_ex, pts_l, boots_syn, popt, data_fs, multinom=True, eps=eps)
    fi.write('Estimated 95% uncerts (with step size '+str(eps)+'): {0}\n'.format(1.96*uncerts_adj[:-1]))
    fi.write('Lower bounds of 95% confidence interval : {0}\n'.format(popt-1.96*uncerts_adj[:-1]))
    fi.write('Upper bounds of 95% confidence interval : {0}\n\n'.format(popt+1.96*uncerts_adj[:-1]))
fi.close()
