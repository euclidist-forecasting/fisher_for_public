""" plfmat.py
    Plot fisher matrices or errors from resutls of the code comparison.
    Markovic, 18/03/2016
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse, traceback, os
import numpy as np
import fishlib as l
import densplotfisher as dp
import itertools as its
import userinterface as ui
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

use_cosmicfish_pylib = True

# import cosmicfish pylib:
if use_cosmicfish_pylib:
    import sys
    cosmicfish_pylib_path = './PlottingScripts/cosmicfish-pyplots'
    sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))
    import cosmicfish_pylib.fisher_matrix as fm
    import cosmicfish_pylib.utilities     as fu
    import cosmicfish_pylib.fisher_operations as fo

# Set all the plot fonts here
LW = 2
FS = 24
matplotlib.rcParams.update({'font.size': FS})
#matplotlib.rcParams['text.usetex'] = True
#Set transparency for lines here
alphatrans=0.6

def processfish(filename, pars, xpars, marg=True, outfilename=''):
    """ Remove and marginalize parameters from Fisher, export Fisher and export paramnames for cosmicfish plots"""
    use_cosmicfish_pylib=True
    data, parlist = l.readfish(filename)
    
    texpars = [l.TeXLabs[par] for par in parlist]
    fidus = [l.FIDUS_BASELINE[par] for par in parlist] 
    fisher = fm.fisher_matrix( fisher_matrix=data, param_names=parlist, param_names_latex=texpars, fiducial=fidus)
    
    if xpars != ['']:
        staypars = [sp for sp in parlist if sp not in xpars]
        fixfish = fo.reshuffle(fisher, staypars)
    
    margfish = fo.marginalise(fixfish, pars)
    if outfilename=='':
        outfilename='_processed'
    outfilename = os.path.splitext(filename)[0]+outfilename+os.path.splitext(filename)[1]
    print("Exporting new Fisher matrix to:  "+outfilename)
    margfish.save_to_file(os.path.splitext(outfilename)[0], simple_header=True)
    return 0

def geterrs(filename, pars, nz, marg=False, notfound=False):
    """ Calculate the marginalised and unmarginalised errors from Fisher matrix file. """

    # Read in the Fisher matrix
    data, parlist = l.readfish(filename)

    # Get the right parameter order dropping any that have been left out
    outpars = l.selectpars(l.zdpars(pars), parlist)
    #print("parlist")
    #print(parlist)
    #print("zpars")
    #print(l.zdpars(pars))
    #print("outpars")
    #print(outpars)
    # initialize cosmicfish fisher matrix:
    if use_cosmicfish_pylib:
        fisher = fm.fisher_matrix( fisher_matrix=data, param_names=parlist )
        if marg:
            data = fo.marginalise( fisher, outpars )
        else:
            data = fo.reshuffle( fisher, outpars )
        eu = np.sqrt(1.0/np.diag(data.fisher_matrix))
        em = np.sqrt(np.diag(data.fisher_matrix_inv))
    else:
        if marg:
            data = np.linalg.inv(l.mreorder(parlist, np.linalg.inv(data), outpars))
        else:
            data = l.mreorder(parlist, data, outpars)
        # Get the marginalised and unmarginalised errors
        eu = np.sqrt(1.0/np.diag(data))
        if np.linalg.det(data) > 0:
            em = np.sqrt(np.diag(np.linalg.inv(data)))
        else:
            print "Skipping " + filename + " marginalised, because Fisher Matrix is singular."
            em = eu

    if notfound:
        return eu, em, outpars, l.selectpars(l.zdpars(pars), parlist, True)[1]
    else:
        return eu, em, outpars

def ploterrs(files, apars, case=4, ts='', marg=True, show=False, zs=l.ZS, alias=None, outpath='./',
             plum=True, yrang=None, cols=l.COLOURS):
    """ Plot the errors vs redshift. """

    # Generic part of name for plot files (if requested)
    if not show: plotfile = os.path.join(outpath, 'errors-comparison-'+str(args.note))

    # Cycle through files and get the errors and the present parameters
    nf = len(files)
    nz = len(zs)

    # Option of legend columns
    if alias is None:
        nc = 1
    elif nf < 6:
        nc = nf
    elif nf >= 6:
        nc = nf/2


    em = nf*[[]]
    eu = nf*[[]]
    dropped = None
    for f,i in zip(files,range(nf)):
        try:
            eu[i], em[i], pars, tmp = geterrs(f, apars, nz, marg, True)
        except Exception as e:
            print traceback.print_exc()
            raise Exception('Error from file ' + f + '.')
        if (dropped is not None) and (len(set(tmp)-set(dropped))>0 or len(set(dropped)-set(tmp))>0):
            raise Exception('Different files contain different sets of parameters!')
        dropped = tmp
    nt = len(pars)

    # Plot differences, not absolute values  np.mean default  np.abs
    mstat_u = np.median(np.array(eu),0)
    mstat_m = np.median(np.array(em),0)
    eurel = 100.0 * (eu-mstat_u)/mstat_u
    emrel = 100.0 * (em-mstat_m)/mstat_m
    max_u = np.max(eurel,0); max_m = np.max(emrel,0)
    min_u = np.min(eurel,0); min_m = np.min(emrel,0)

    # Save errors to file if plots are being saved too
    if not show:
        # Save the errors in a file
        header = ' '.join(pars) + '\n'
        header += 'The columns are in the order of the above parameters.\n'
        header += 'The rows are in groups of four: 1. unmarginalised errors, 2.marginalised errors, 3. marginalised relative error differences w.r.t median, 4. unmarginalised relative error differences w.r.t median, for each of the below files:'
        for f in files:
            header += '\n' + f
        errfile = os.path.join(outpath, str(args.note)+'-errors.dat')
        allerr = []
        for i in range(nf):
            allerr.append(eu[i])
            allerr.append(em[i])
            allerr.append(emrel[i])
            allerr.append(eurel[i])
        #allerr = np.array(np.transpose(allerr))
        allerr = np.array(allerr)
        #print allerr.shape
        np.savetxt(errfile, allerr, fmt='%.10e', header=header)
        print 'Errors saved to file ' + errfile + '.'


    ###
    # First make plots for all the Z-DEPENDENT parameters against redshift
    #    save the locations of the shape parameters for later
    ###
    iz_start = 0
    ishape = []
    shape_pars = []
    for par in apars:

        # make a brand new plot for each z-dependent parameter
        if par in l.ZPARS:
            fig = plt.figure(figsize=(10,6))
            ax = fig.add_subplot(111)
            ax.set_position([0.1,0.1,0.8,0.8])
            labpar=l.TeXLabs[par]
        elif par not in pars:
            continue
        elif par in l.PARS:
            # Save for later - go down the list by 1 row only
            ishape.append(iz_start)
            shape_pars.append(l.TeXLabs[par])
            iz_start+=1
            continue

        for i in range(nf):

            # Set alias in the legend if given
            if alias is None:
                lbl = l.shorten(args.files[i].split('/')[-1],40)
            else:
                lbl = alias[i]


            # Plot the marginalised errors
            ax.plot(zs, emrel[i][iz_start:iz_start+nz], ls=':', marker=".",
                    ms=LW*12, c=cols[i], lw=LW, label=lbl, alpha=alphatrans)

            # Plot the unmarginalized errors
            if plum==True:
                ax.plot(zs, eurel[i][iz_start:iz_start+nz], ls='-', ms=LW*6, mew=6, c=cols[i], lw=LW, alpha=alphatrans)

        plt.fill_between(zs, min_m[iz_start:iz_start+nz], max_m[iz_start:iz_start+nz], interpolate=True, facecolor='0.9', edgecolor='0.9', alpha=0.5, linewidth=0.0)
        plt.fill_between(zs, min_u[iz_start:iz_start+nz], max_u[iz_start:iz_start+nz], interpolate=True, facecolor='0.3', edgecolor='0.3', alpha=0.5, linewidth=0.0)
        iz_start+=nz

        darkgreypatch = mpatches.Patch(color='0.3', alpha=0.5)
        # Prettyfy the z-dependent plots and save if requested:
        #
        dashart = plt.Line2D((0,1),(0,0), color='0.3', lw=1.5, linestyle='-', ms=LW*6, mew=6, alpha=alphatrans)
        solidart = plt.Line2D((0,1),(0,0), color='k', lw=1.5, linestyle=':', marker=".", markerfacecolor='0.9', ms=LW*12, alpha=alphatrans)
        if plum==True:
            leg1 = ax.legend([solidart, darkgreypatch], ['marg.', 'unmarg.'], loc='best', ncol=2, prop={'size':14}, handlelength=2, numpoints=1)
        else:
            dashart = plt.Line2D((0,1),(0,0), color='0.3', marker="s", ms=LW*6, lw=0, markeredgewidth=0.0, alpha=alphatrans)
            leg1 = ax.legend([dashart, solidart], ['unmarg.', 'marg.'], loc='best', ncol=2, prop={'size':14}, handlelength=2, numpoints=1)
        #
        #ax.axhline(y=0.0, ls=':', c='k')
        ax.add_artist(leg1)
        plt.legend(bbox_to_anchor=(1, 1), loc='lower right', ncol=nc, prop={'size':12})
        plt.ylabel(r'% differences on ' +r'$\sigma_{'+labpar+r'}$')
        plt.xlabel('z')
        plt.xlim([min(zs)-0.05, max(zs)+0.05])
        ax.tick_params(axis='x', pad=10)
        
        majorLocator = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter('%.1f')
        minorLocator = MultipleLocator(0.05) 
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)
        
        majorLocator = MultipleLocator(2.0)
        minorLocator = MultipleLocator(1.0) 
        majorFormatter = FormatStrFormatter('%.0f')

        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        if yrang is None:
            yymin, yymax = plt.ylim()
            yymax = np.max(np.abs([yymin,yymax]))
            if yymax < 0.05:
                yymax = 0.05
            if yymax < 0.1:
                yymax = 0.1
            elif yymax < 1.0:
                yymax = 1.0
            elif yymax < 3.0:
                yymax = 3.0
            elif yymax < 6.0:
                yymax = 6.0
            elif yymax < 10.0:
                yymax = 10.0
            else:
                yymax = 1.05*yymax
            plt.ylim([-yymax, yymax])
        else:
            plt.ylim([-yrang, yrang])
        if not show: fig.savefig(plotfile+'-'+par,dpi=400,bbox_inches='tight')

    ###
    # Now plot the redshift-independent (SHAPE) parameters if there are any
    ###
    if len(ishape)>0:
        fig = plt.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        ax.set_position([0.1,0.1,0.8,0.8])

        for i in range(nf):

            # Set alias in the legend if given
            if alias is None:
                lbl = l.shorten(args.files[i].split('/')[-1],40)
            else:
                lbl = alias[i]


            # Plot the marginalised errors
            ax.plot(ishape, emrel[i][ishape], '.', c=cols[i], ms=LW*12, alpha=alphatrans, label=lbl)

            # Plot the unmarginalized errors
            if plum==True:
                ax.plot(ishape, eurel[i][ishape], linestyle='', c=cols[i], marker='', ms=0, mew=6, alpha=alphatrans)

        # Plot the shaded regions
        plt.fill_between(ishape, min_m[ishape], max_m[ishape], interpolate=True, facecolor='0.9', edgecolor='0.9', alpha=0.5, linewidth=0.0)
        plt.fill_between(ishape, min_u[ishape], max_u[ishape], interpolate=True, facecolor='0.3', edgecolor='0.3', alpha=0.5, linewidth=0.0)

        # Prettyfy the shape plot and save if requested:
        #
        #plt.figure(fig.number)
        pointart = plt.Line2D((0,1),(0,0), color='k', marker='.', lw=1, ms=LW*12, linestyle='None')
        crossart = plt.Line2D((0,1),(0,0), color='k', marker='x', lw=1, ms=LW*6, mew=6, linestyle='None')
        darkgreypatch = mpatches.Patch(color='0.3', alpha=0.5)
        lightgreypatch = mpatches.Patch(color='0.9', alpha=0.5)
        if plum==True: # Which code legend
            leg2=ax.legend([lightgreypatch, darkgreypatch], ['marg.', 'unmarg.'], loc='upper right', ncol=2, prop={'size':14})
        else:
            dashart = plt.Line2D((0,1),(0,0), color='0.3', marker="s", ms=LW*6, lw=0, markeredgewidth=0.0, alpha=alphatrans)
            leg2 = ax.legend([dashart, solidart], ['unmarg.', 'marg.'], loc='best', ncol=2, prop={'size':14}, handlelength=1, numpoints=1)
        #
        plt.xticks(ishape, shape_pars)
        ax.legend(bbox_to_anchor=(1, 1.05), loc='lower right', ncol=nc, prop={'size':14}, handlelength=2, numpoints=1)
        ax.add_artist(leg2)
        #ax.axhline(y=0.0, ls=':', c='k', alpha=0.2) # we don't want the zero-line
        plt.ylabel(r'% differences on ' +r'$\sigma_i$')
        plt.xlim([min(ishape)-0.2, max(ishape)+0.2])
        ax.tick_params(axis='x', pad=10)
        #
        # Set the y-range
        if yrang is None:
            yymin, yymax = plt.ylim()
            yymax = np.max(np.abs([yymin,yymax]))
            if yymax < 0.3:
                yymax = 0.3
                locat = 0.05
            elif yymax < 1.0:
                yymax = 1.0
                locat = 0.1
            elif yymax < 3.0:
                yymax = 3.0
                locat = 1.0
            #elif yymax < 6.0:   #remove this bound, so all comparison plots of paper are from -10 to +10
            #    yymax = 6.0
            #    locat = 2.0
            elif yymax < 10.0:
                yymax = 10.0
                locat = 5.0
            else:
                yymax = 1.05*yymax
                locat = yymax//5
            plt.ylim([-yymax, yymax])
        else:
            plt.ylim([-yrang, yrang])
        #
        majorLocator = MultipleLocator(locat)
        minorLocator = MultipleLocator(locat/2) 
        majorFormatter = FormatStrFormatter('%.0f')

        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        if not show: fig.savefig(plotfile,dpi=400,bbox_inches='tight')

    return 0

def plotfmat(files, pars, case=4, ts='', marg=True, show=False, outpath='./'):
    """ 2D plot of the Fisher matrix. """

    plotfile = 'fisher-matrix-'+ts+'-case'+str(case)+'-note-'+str(args.note)+'-file'
    plotpars = l.selectpars(l.PARS+l.parszd(), pars)

    nf = len(files)
    nt = len(plotpars)
    nz = len(l.ZS)
    matchings=None
    matchind = None
    if nf > 1:
        matchings = list(its.combinations(files, 2)) #produces all possible combinations of files
        matchind = list(its.combinations(range(nf), 2))
    elif nf==1:
        matchings = [(files[0], files[0])]
        matchind = [(i,i) for i in range(1)]
    if marg==True:
        margstr = 'marginalised'
    else:
        margstr = 'unmarginalised'
    # Cycle through files
    for f12, ij in zip(matchings, matchind):

        # Read in the Fisher matrix & process it
        f1 = f12[0]
        f2 = f12[1]
        indid = '-'.join(str(k) for k in ij)

        data1, parlist1 = l.readfish(f1)
        data2, parlist2 = l.readfish(f2)

        ftitle1 = f1.split('/')[-1]
        ftitle2 = f2.split('/')[-1]
        print "   time stamp:\n\t" + ts
        if ftitle1 == ftitle2:
            plottit = ftitle1
            print "   file:\n\t" + ftitle1
        else:
            plottit = ftitle1+' \n '+ftitle2
            print "   files:\n\t" + ftitle1 + "\n\t" + ftitle2
        if len(data1)==nt:
            plottit += ' \nall parameters'
        else:
            plottit += ' \n'+margstr+' over missing parameters'

        if not marg:
            data1 = l.mreorder(parlist1, data1, plotpars)
            data2 = l.mreorder(parlist2, data2, plotpars)
        else:
            data1 = l.marginalise(parlist1, data1, keep=plotpars)
            data2 = l.marginalise(parlist2, data2, keep=plotpars)
        if nt > 6:
            printperc = False  #if matrix to be plotted is bigger than 6x6, then don't plot numbers on it
            fs = 4
        else:
            printperc = True
            fs = 6
        # Plot the whole thing
        #plt.imshow(data, interpolation='nearest')
        #plt.imshow(err_unm, interpolation='nearest')
        #if not show: fig.savefig(plotfile+'-'+str(i),dpi=400,bbox_inches='tight')

        dp.print_matrix_with_values(data1, data2,
                                    out_path=outpath,
                                    pars=plotpars,
                                    exportfile=(not show),
                                    png_name=plotfile+'-'+margstr+'-comb_'+str(indid),
                                    plot_title=plottit,
                                    submatrix_case="full",
                                    print_percentages=printperc,
                                    font_size=fs )

    return 0

def check_fmat(args, nz=len(l.ZS)):
    """ Run a quick check that the input arguments and the file given are consistent."""

    # First check if given files or folder(s)
    addf = []
    addal = []
    addcol = []
    for fpath in args.files:
        if os.path.isdir(fpath):
            addf.extend([os.path.join(fpath, f) for f in os.listdir(fpath) if (os.path.isfile(os.path.join(fpath, f)) and ('.dat' in f or '.txt' in f))])
            args.files.remove(fpath)
    args.files.extend(addf)   #appends to present file list all files within folder that match .dat or .txt
    #print(args.files)
    #print(args.alias)
    #print(args.colours)
    for fpath, falias, fcol in zip(args.files, args.alias, args.colours):
        if os.path.isfile(fpath):
            addf.append(fpath)
            addal.append(falias)
            addcol.append(fcol)
        else:
            print('warning: '+os.path.basename(fpath)+' not found in folder!')
    if len(addf)==0:
        return "No matching filename found in folder"
    args.files = addf   #overwrites args.files with only valid file names, and matches the aliases and colours
    args.alias = addal
    args.colours = addcol

    # Cycle through files
    for f in args.files:
        data, parlist = l.readfish(f)
        if data.shape[0]!=data.shape[1]:
            return 'Matrix in ' + f + ' is not square!'
        if args.case < 3 and data.shape[0] > 56:
            return "You have included P_shot in the Fisher matrix even though its set to 0!"
        if args.case < 5 and data.shape[0] % nz > 0:
            return 'Either you are using a wrong number of redshift bins, the wrong case (c=' +\
                   str(args.case) + '), or you have redshift independent parameters in your matrix.'
        if args.case == 5:
            l.PARS = l.shapePARS
        if args.case == 6:
            l.PARS = l.flatPARS
        if args.case == 7:
            l.PARS = l.nonflatPARS
        if args.case == 8:
            l.PARS = l.flatgammaPARS
        if args.case == 9:
            l.PARS = l.nonflatgammaPARS
        if args.case == 10:
            l.PARS = l.shapePARS+l.NLPARS
        if args.case == 60:
            l.PARS = l.flatPARS+l.WLNUISANCE
        if args.case == 70:
            l.PARS = l.nonflatPARS+l.WLNUISANCE
        if args.case == 80:
            l.PARS = l.flatgammaPARS+l.WLNUISANCE
        if args.case == 90:
            l.PARS = l.nonflatgammaPARS+l.WLNUISANCE
        if data.shape[0] < len(args.pars):
            return 'You are asking for more parameters than how many are in the matrix.'
    return None


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Make Fisher Matrix plots for the IST SG3 (Galaxy Clustering) code comparison.")
    parser.add_argument('files', type=str, nargs='+', help='a list of files to plot')
    parser.add_argument('-o', '--outpath', type=str, default='./', help='folder where to save the plots to')
    parser.add_argument('-a', '--alias', type=str, nargs='+', help='aliases for the legend (order should corespond to file order)')
    parser.add_argument('-r', '--colours', type=str, nargs='+', help='colours for the legend (order should corespond to file order)')
    parser.add_argument('-c', '--case', type=int, default=5, choices=l.CASES, help='which case you are inputting')
    parser.add_argument("-s", "--show", action = 'store_true', default=False, help='plotting flag - set true to show instead of saving')
    parser.add_argument("--notstamp_outdir", action = 'store_false', default=True, help='Flag, if not present, create a subfolder in outpath with the timestamp name to store the output')
    parser.add_argument("-g", "--plmargonly", action = 'store_true', default=False, help='Flag: If present, plot only the marginalised errors for the parameters.')
    parser.add_argument("-m", "--marg", action = 'store_true', default=False, help='marginalise flag for plotting the matrix. does noting if error flag is on')
    parser.add_argument("-p", "--pars", default=l.PARS+l.NLPARS+l.ZPARS, nargs='+', help='list of parameters to plot - must go after the filename')
    parser.add_argument("-x", "--fixpars", default=[''], nargs='+', help='list of parameters to remove from fisher')
    parser.add_argument("-y", "--yrange", type=float, default=None, help='(default: set by maximum value in each plot) half of y-range in percentages')
    parser.add_argument("-n", "--note", default='', help='note for the log file and for the filename of plots and folders')
    pars = parser.add_mutually_exclusive_group(required=False)
    pars.add_argument("-e", "--errors", action = 'store_true', default=False, help="(default: True) flag specifying that files contain the full Fisher matrix created for the corresponding case")
    pars.add_argument("-f", "--fmat", action = 'store_true', default=False, help="flag specifying that files contain the full Fisher matrix and that the whole matrix should be plotted for each")
    pars.add_argument("-i", "--imat", action = 'store_true', default=False, help="flag specifying that the fisher matrix and its bounds are plotted")
    args = parser.parse_args()
    if not args.errors and not args.fmat and not args.imat: args.errors = True

    ts = l.get_timestamp()
    if args.notstamp_outdir == True:
        ts = ts+'-'+args.note    
        args.outpath = os.path.join(args.outpath, ts)
        print args.outpath


    if args.alias is None or len(args.alias) != len(args.files):
        args.alias = [l.shorten(args.files[ii].split('/')[-1],40) for ii in range(len(args.files))]
    if args.colours is None or len(args.colours) < len(args.files):
        args.colours = [l.COLOURS[ii] for ii in range(len(args.files))]
    # First run a quick check of the input, throw exception if it looks wrong
    err = check_fmat(args)
    if err is not None:
        raise Exception(err)

    # after checks passed, creates directory if non-existent and if show option is not given
    if not args.show:
        ui.mkdirp(args.outpath)
    
    if args.plmargonly==True:
        plotunmarg=False
        plotstr="Plotting marginalised errors."
    else:
        plotunmarg=True
        plotstr="Plotting marginalised and unmarginalised errors."
    # Create the plots
    if args.fmat:
        print "Plotting Fisher matrices."
        O = plotfmat(args.files, args.pars, args.case, ts, args.marg, args.show, args.outpath)
    elif args.errors:
        print plotstr
        O = ploterrs(args.files, args.pars, args.case, ts, args.marg, args.show, l.ZS, args.alias,
                     args.outpath, plotunmarg, args.yrange, args.colours)
    elif args.imat:
        for fi in args.files:
            O = processfish(fi, args.pars, args.fixpars, outfilename=args.note)
    # Save & log or show
    if args.show:
        plt.show()
    else:
        logfile = os.path.join(args.outpath, str(args.note)+'.log')
        l.savelog(args, logfile, info=args.note)
