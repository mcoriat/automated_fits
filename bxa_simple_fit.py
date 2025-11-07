import xspec
from xspec import *
import bxa.xspec as bxa
import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.io import fits

# Load spectrum and background
spectrum = Spectrum("3067718060100029_pn.grp")
spectrum.background = "P0760940101PNS003BGSPEC0017.FTZ"

# Load RMF and ARF files
AllData(1).response = "epn_e3_ef20_sdY6.rmf"
AllData(1).multiresponse[0].arf = "P0760940101PNS003SRCARF0017.FTZ"

# Ignore channels outside desired energy range (0.3–10 keV)
AllData.ignore("**-0.3 10.0-**")

# Define XSPEC model (absorbed power-law, fixed redshift)
model = Model("phabs*zpowerlw")
model.zpowerlw.Redshift = 1.0  # fixed redshift

# Set parameter ranges explicitly on XSPEC parameters
model.phabs.nH.values = "0.05,,0.001,0.001,10.0,10.0"
model.zpowerlw.PhoIndex.values = "2.0,,1.0,1.0,3.0,3.0"
model.zpowerlw.norm.values = "1e-4,,1e-6,1e-6,1e-2,1e-2"

# Explicit transformations (priors) in oldest BXA syntax (corrected!)
transformations = [
    bxa.create_uniform_prior_for(model, model.phabs.nH),
    bxa.create_uniform_prior_for(model, model.zpowerlw.PhoIndex),
    bxa.create_loguniform_prior_for(model, model.zpowerlw.norm)
]

# Set the required Cash statistic explicitly
Fit.statMethod = "cstat"

# Run BXA solver
solver = bxa.BXASolver(transformations=transformations, outputfiles_basename='bxa_fit_result/')
results = solver.run(resume=False)

# Plot posterior corner plot
#solver.plot()


# Print final best-fit parameters
#print("Best-fit parameters (posterior median):")
#for par, median in zip(solver.paramnames, results["posterior_median"]):
#    print(f"{par}: {median:.4e}")




# Load posterior samples from the BXA output directory
with fits.open("bxa_fit_result/chain.fits") as hdul:
    samples = hdul[1].data
    samples_array = np.column_stack([samples[name] for name in samples.names])

# Get sample array and labels
labels = solver.paramnames
stds = np.std(samples_array, axis=0)
valid_cols = stds > 0  # only include parameters with non-zero spread
filtered_samples = samples_array[:, valid_cols]
filtered_labels = [label for i, label in enumerate(labels) if valid_cols[i]]

# Plot only parameters with dynamic range
fig = corner.corner(filtered_samples, labels=filtered_labels, show_titles=True, title_fmt=".3e")
fig.savefig("bxa_fit_result/corner_plot.png")
plt.show()




# Load posterior samples
#samples = np.loadtxt("bxa_fit_result/posterior.txt")

# Compute medians
medians = np.median(samples_array, axis=0)
print("Best-fit parameters (posterior median):")
for i, par in enumerate(labels):
    print(f"{par}: {medians[i]:.4e}")

