""" densplotfisher.py
Functions to plot entries of a Fisher matrix or the percentage difference between two Fisher matrices
Santiago Casas 20/08/2016
"""


import numpy as np
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.colors import Normalize
import os.path

# this function selects a submatrix of the Fisher matrix, according to a string keyword.
def submatr_cases(stri, st):
    sfull=slice(None)
    sshap=slice(0,st)
    szdep=slice(st,None)
    return {
        "full"      :(sfull,sfull),
        "shape"     :(sshap,sshap),
        "zobs"      :(szdep,szdep),
        "shapezcorr":(sshap,szdep)
    }.get(stri,(sfull,sfull))

# This is to get white at zero
# Copied from here: http://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


#this function prints a colored MatrixPlot comparing two Fisher matrices. If chosen, it prints on each square the percentage difference.
#If a single Matrix is given, then its entry values are plotted. 
def print_matrix_with_values(matrix1, matrix2, 
                             out_path='./', 
                             png_name='matrix.png', 
                             plot_title='',
                             font_size=4, 
                             font_size_title=8, 
                             cmap='gray', 
                             print_percentages=True, 
                             slicestart=4, 
                             submatrix_case="full", 
                             exportfile=True, 
                             pars=None,
                             get_submat=False):
    
    # Pretty fonts
    font_ticks = {'size': font_size}
    matplotlib.rc('font', **font_ticks)

    # Create new plot
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.set_title(plot_title, fontsize=font_size_title, va='bottom', )
    
    # Extract part of the matrix if not the full thing is needed
    if submatrix_case != 'full' and get_submat:
        sl = submatr_cases(submatrix_case, slicestart)        
        print 'slicing matrix:\n\t' + str(sl)
        matrix1 = matrix1[sl]
        matrix2 = matrix2[sl]

    # Assign the matrix to plot: either fishermatrix or the comparison between two of them
    dsign = np.multiply(np.sign(matrix1),np.sign(matrix2))
    if np.array_equal(matrix1,matrix2):
        # Result assignment and text label creation
        matrix=matrix1
        formatstring = '{:-.2e}'
        fontnum = font_size
        addtxt = ''
    else:
        # Calculate the comparison metric and create text labels
        rho1 = 0.0*np.array(matrix1)
        rho2 = 0.0*np.array(matrix2)
        for i,row in enumerate(rho1):
            rho1[:,i] = np.divide(np.ravel(matrix1[:,i]), \
                        np.sqrt(np.diag(matrix1)[i] * np.diag(matrix1)))
            rho2[:,i] = np.divide(np.ravel(matrix2[:,i]), \
                        np.sqrt(np.diag(matrix2)[i] * np.diag(matrix2)))
        matrix = np.divide(rho1+1e-16, rho2+1e-16)
        formatstring = '{:.2f}'
        fontnum = font_size
        addtxt = 'r1/r2='

    # Plot color matrix:
    cax = ax.matshow(matrix, cmap=pl.cm.Spectral, norm=MidpointNormalize(midpoint=1.0, vmax=3.0, vmin=-1.0))
    pl.xticks(range(matrix.shape[0]), pars, rotation=80)
    pl.yticks(range(matrix.shape[1]), pars)
    cbar = fig.colorbar(cax, extend='max')
    ax.set_aspect('equal')
    cbar.ax.tick_params(labelsize=10)
    if addtxt: cbar.ax.set_ylabel('ratio of correlation coefficients', rotation=270, fontsize=12, labelpad=0.5)
    cbar.ax.yaxis.set_label_position('left') 

    # Text portion:
    if print_percentages==True:
        ind_x, ind_y = np.arange(matrix.shape[0]), np.arange(matrix.shape[1])
        x, y = np.meshgrid(ind_x, ind_y)
        for x_val, y_val in zip(x.flatten(), y.flatten()):
            txt = addtxt+formatstring.format(matrix[x_val,y_val])
            if dsign[x_val,y_val]< 0: txt += ' (-!)'
            ax.text(x_val, y_val, txt, va='center', ha='center', fontsize=fontnum, color="black", bbox=dict(facecolor='white', edgecolor='black', boxstyle='round' ))
    if exportfile==True:
        output_png_path = os.path.join(out_path,submatrix_case.title()+'_'+png_name)
        fig.savefig(output_png_path, dpi=400, bbox_inches='tight')
    else:
        pl.show()
        
# Test the code
if __name__ == '__main__':

    import fishlib as l

    # What pars to plot
    sub_case = "shape"
    #plotpars = ['wm','wb','ns','h']
    plotpars = ['lnH']
    plotpars = l.selectpars(l.PARS+l.parszd(), plotpars)

    # Get the files
    PATH_IN = "./"    
    file_name2 = "alkistis-case5-with-files-v0.2.txt"
    file_name1 = "sapone_shape_zobs_ist_case5_pk.dat"

    # Get the FM data
    matrix1, pars1 = l.readfish(PATH_IN+file_name1)
    matrix1 = l.mreorder(pars1, matrix1, plotpars)
    matrix2, pars2 = l.readfish(PATH_IN+file_name2)
    matrix2 = l.mreorder(pars2, matrix2, plotpars)
    sub_case='full' # this is a hack, because we are already selecting the submatrix here
    
    # Title for plot
    plottit = "Alkistis vs. Santiago input"
    
    # Make plot        
    print_matrix_with_values(matrix1, matrix2, out_path=PATH_IN, exportfile=False, 
                             plot_title=plottit, font_size='6', font_size_title='10', 
                             print_percentages=False, submatrix_case=sub_case,
                             pars=plotpars)
