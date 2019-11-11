# Fisher_for_public

*All the Sections or Figures mentioned in this readme refer to the IST:F paper https://arxiv.org/abs/1910.09273*

**CONTENT OF THIS REPOSITORY**

This repository contains representative Fisher matrices obtained following the recipe of the IST:F forecast paper.


One Euclid IST:F reference Fisher matrix is provided for each probe and cosmology.
These matrices are in the `All_results` folder, divided in the `optimistic` and `pessimistic` subfolders, corresponding to the 2
sets of specifications used in this paper. These subfolders, further contain a `flat` and `non-flat` folder each. Within these, the files
containing the Fisher matrices follow the naming convention

`EuclidISTF_Observables_CosmoModel_SpecificationsCase.txt`

For example, a Fisher matrix whose file name is

`EuclidISTF_GCph_WL_XC_w0wa_flat_optimistic.txt`

corresponds to the combined Fisher matrix for photometric galaxy clustering, weak lensing and their cross correlation, obtained in
the flat case for the flat w0,wa cosmology, using the optimistic set of specifications. Each Fisher matrix contains a header listing the
cosmological parameter of the matrix and their order.

In addition, the repository contains a set of python scripts that can be used to plot the contour ellipses and compare them with
user-provided matrices.

**HOW TO PRODUCE A FISHER MATRIX USING IST:F INPUT**

The Euclid matrices contained in the public repository are obtained following the recipes described in section 3 for the cases
analysed in section 4 and section 5.

Users who wish to validate their own code can build their Fisher matrices following these recipes, using the input files for
cosmological quantities that are also included within the repository. These input files have been generated using CAMB, version
August 2018.

The input files containing the cosmological quantities used to create the Fisher matrices can be found at 
http://pc.cd/0wJotalK .
The link contains two folders: GC and WL.


*GC folder:*

the `GC` folder contains the matter power spectrum P(k,z) and growth function D(z) generated using CAMB. These
quantities are available in the zip file `ISTF_GC.zip`. The parent folder contains the growth function D(z) and the growth rate f(z).
The subfolder `Pk` contains the Pm(k,z), while the subfolder Pk-nw contains the Pnw(k,z) that appears in Equation 80.

Both these folders are organized as follows:

- in the fiducial folder, the P(k,z) for the fiducial cosmology are contained, with one file for each redshift at which this spectrum
is computed. Each file has a suffix `z_ii`, where `ii` runs from 00 to 03 corresponding to the 4 redshift bins used for the Euclid
GCsp probe `z = {1.0,1.2,1.4,1.65}`.
- in the `par_step_eps_1p0E-2` folders, with `par` being the name of the parameter and `step` indicating a positive (pl) or
negative (mn) step, the power spectra for each step in cosmological parameters are contained, following the same structure of
the fiducial one (the steps amplitudes are relative; for each parameter the absolute step will be `fidpar(1 + eps)` with `fidpar`
the fiducial value of the parameter and `eps` the value of the relative step).
All the power spectrum files have the column structure

<pre>
k [h=Mpc]      linear P(k,z)/sigma8^2(z) [Mpc^3/h^3]      sigma8^2(z)
</pre>


while the f(z) and D(z) files have the structure

<pre>
z       f(z)        D(z)
</pre>

*WL folder:* 

the `WL` folder contains the matter power spectrum P(k,z), generated using CAMB for the fiducial model and for
each step in cosmological parameters used to compute numerical derivatives, together with a `scaledmeanlum-E2Sa.dat` file,
which contains `z` (first column) and the ratio `<L>(z)/L*(z)` (second column) described in Equation 109.

The matter power spectra are available in two zip files:

- `ISTF_WL_Flat.zip`
- `ISTF_WL_NonFlat.zip`

The first contains the matter power spectrum P(k,z) for the flat case, while the second contains the same for the non-flat case.
In each of these two zip files, the fiducial power spectrum can be found in the file `pkz-Fiducial.txt`, while the spectra relative
to the steps in cosmological parameters are contained in subfolders whose name describe the step taken, e.g. `ns_pl_eps_5p0E-2`
indicates a step on ns in the positive direction of amplitude `5*10^-2` (the steps amplitudes are relative; for each parameter the
absolute step will be `fidpar(1 + eps)` with fidpar the fiducial value of the parameter and `eps` the value of the relative step).

Notice that with respect to the GC input files, here many more steps for numerical derivatives are computed; this is because most
of the IST:F codes for WL use the derivative scheme described in subsection 4.3.

All the P(k,z) files have the column structure

<pre>
z     k [h=Mpc]     P(k,z) linear [Mpc^3/h^3]     P(k,z) nonlinear [Mpc^3/h^3]
</pre>


**HOW TO INCLUDE AN EXTERNAL MATRIX IN THE REPOSITORY**

In order to compare their own Fisher code with the IST:F results, users will need to reproduce the cases described here and add their
matrices in the corresponding folder. The filenames of the new matrices, must follow the convention of the repository, i.e. with a
prefix stating the name of the code used, and a suffix containing the observables included in the matrix, the cosmological model used
and the specifications case. The full description of the cases, with all the corresponding specifications used to obtain the matrices,
can be found in section 4.

Furthermore, each matrix file provided by users must contain the header listing the cosmological parameters and their order.
The matrices provided by the user must be placed in the subfolders of the corresponding cases.

**HOW TO COMPARE AN EXTERNAL FISHER MATRIX**

The python script performing the IST:F comparison can be launched from the parent folder using

`python comparison.py [-o OBSERVABLES ] [-f USERFISHER ]`

The `-o` and `-f` arguments are optional.

`USERFISHER` are the labels of the user provided Fisher matrices. Notice that the user's matrices must be placed in the correct
folders of the repository and must follow the filename notation described above, with `USERFISHER` the prefix of the filename, i.e.
the code name. The script will produce plots, similar to the comparison plots of section 4, only for those cases for which an external
matrix to compare with is found. In case no `USERFISHER` argument is provided, a test matrix will be used. This test matrix is
contained within the repository only for the Cross Correlation comparison cases (see subsection 4.4).

The `OBSERVABLES` can be one or more of the cases considered in section 4, i.e.

- GC for spectrocopic Galaxy Clustering
- WL for Weak Lensing
- XC for Probe Combination and Cross Correlation

The script can perform one, two or all cases at once. The latter is the default selection if no observable is specified. Results will
be placed in results folders in each of the subcases locations, e.g.

`All_results/optimistic/flat/results_XC`

The comparison tool consists of:

- a common trunk, comparison.py found in the parent folder.
- 3 apps (appGC.py, appWL.py, appXC.py) found in PlottingScripts/ handling the different comparisons
- a set of plotting and analysis libraries found in PlottingScripts/

As an example, let's assume that a user wants to compare a WL matrix computed in the optimistic case, for a flat `w_0,w_a` cosmology. The user will add the matrix as

`All_results/optimistic/flat/MyName_WL_w0wa_flat_optimistic.txt`

Running the script with the command

`python comparison.py -o WL -f MyName`

will compare the new matrix with the Euclid IST:F for this specific case and cosmology, which will be found in

`All_results/optimistic/flat/results_WL` 

Notice that, as said above, the script does not produce comparison results for those cases in which no external matrix is provided. In the example above, where only the WL, optimistic matrix for `w_0,w_a` is provided, the script will produce the results only for this specific case also running

`python comparison.py -f MyName`

**HOW TO OBTAIN RESULTS FOR THE IST:F AND EXTERNAL MATRICES**

The python script performing the IST:F ellipses plot can be launched from the parent folder using


`python ellipses.py [-o OBSERVABLES ] [-f USERFISHER ] [-n USERNAME ] [-c USERCOLOR]`

The `-o`, `-f`, `-n` and `-c` arguments are optional.

In case no arguments are provided,  the script will reproduce Fig.10 from the paper.

With the option `-o` a list of `OBSERVABLES` can be specified, separated by spaces. Available options: 

- `All` produces a plot like Fig.10 for different observables.
- `XC` produces a plot focusing on XC correlations, like Fig.11.
- `GC` produces a plot containing only GCs contours.
- `WL` produces a plot containing only WL contours.

With the option `-f`, a `USERFISHER` can be specified, which is the relative path to the Fisher matrix provided by the user. Notice that the user's matrices must be placed in the correct folders of the repository and must follow the filename notation described above. **Warning: while in the comparison script case the user can provide multiple external Fisher matrices, the results scripts onle accept a single external matrix.**

If an external user Fisher matrix is provided, the script will overplot it on top of the other ellipses using  a color and a label name provided by the user. 
Default color: `black`, default name: `USER`.

The user can change the color and the label name of its Fisher matrix, with the option `-c` and `-n`, respectively.

Ellipses plots will be placed in the following folder, e.g.

`All_results/optimistic/flat/results_ellipses`

Additionally to the elliptical contours, the routine also produces a file called:

`w0waCDM-flat-optimistic_bounds.txt`

containing the relative errors on all the parameters of interest and a file

`w0waCDM-flat-optimistic_FoMs.txt`

containing the FoM for `w_0,w_a` corresponding to each of the Fisher matrices which were plotted.

The `.ini` files inside the directory `PlottingScripts/` store all the parameters of the plot, and provide a way for the user to gain finer control of the plotting routine, if needed. 

