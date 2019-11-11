#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2017 by the CosmicFish authors
#
# The CosmicFish code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CosmicFish distribution.
#
#----------------------------------------------------------------------------------------

"""

Simple Python code to perform analysis of Fisher matrices 

The ouput will be a file with bounds and FoMS


"""

# ***************************************************************************************

__version__ = '1.0' #: version of the application

# ***************************************************************************************

""" Hard coded options """

# ***************************************************************************************

# import first dependencies:
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines    as mlines
import numpy as np
import argparse
import math
import sys
import os
import copy
import itertools as it
import ConfigParser

# get the path of the application and the CosmicFish library:
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

# import the CosmicFish pylib
import cosmicfish_pylib.utilities            as fu
import cosmicfish_pylib.colors               as fc
import cosmicfish_pylib.fisher_matrix        as fm
import cosmicfish_pylib.fisher_derived       as fd
import cosmicfish_pylib.fisher_operations    as fo
import cosmicfish_pylib.fisher_plot_settings as fps
import cosmicfish_pylib.fisher_plot_analysis as fpa
import cosmicfish_pylib.fisher_plot          as fp


# ***************************************************************************************

def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def get_with_default(section,name, default=''):
    try:
        return Config.get(section,name)
    except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
        return default

# protection against importing:
if __name__ == "__main__":
    
    # parse command line arguments:
    parser = argparse.ArgumentParser(description='Analysis tool for plot and bounds')
    # parse file names:
    parser.add_argument('inifile', metavar='inifile', type=str, nargs='+',
                         help='file with a list of instructions')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
 

    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header('Global analysis app version '+__version__)

    # process input arguments:
    inifile          = args.inifile

    #function used to deal with ini sections
    def ConfigSectionMap(section):
        dict1 = {}
        options = Config.options(section)
        for option in options:
            try:
                dict1[option] = Config.get(section, option)
                if dict1[option] == -1:
                   DebugPrint("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1
    
    #initializing and reading the config file
    Config = ConfigParser.ConfigParser()
    Config.read(inifile)

    #Reading general options
    outroot   = ConfigSectionMap("General Options")['outroot']
    files     = Config.get("General Options", "fishers").split("\n")
    derived   = Config.getboolean("General Options", "derived")
    sum_fish1  = get_with_default("General Options", "sum_fish1").split("\n")
    sum_fish2  = get_with_default("General Options", "sum_fish2").split("\n")
    sum_fish3  = get_with_default("General Options", "sum_fish3").split("\n")
    sum_fish4  = get_with_default("General Options", "sum_fish4").split("\n")
    sum_fish5  = get_with_default("General Options", "sum_fish5").split("\n")
    usercolors  = get_with_default("General Options", "colors").split("\n")
    eliminate = Config.getboolean("General Options", "eliminate")
    fishnames = Config.get("General Options", "names").split("\n")
    pform     = '.pdf'
    sum_fish = [sum_fish1, sum_fish2, sum_fish3, sum_fish4, sum_fish5]
    summax = len(sum_fish)
    #General screen output
    if not args.quiet:
       print 'GENERAL OPTIONS:'
       print ' Output root='+outroot
       print ' Using derived parameters='+str(derived)
       print ' Eliminate rather than marginalize='+str(eliminate)
       print ' ---------------------------------'
       print ' Bounds from these matrices will be computed:'
       for elem in files:
           print elem
       for si in range(summax):
           if sum_fish[si][0]:
               print '#1 Also the combination of these will be computed:'
               for elem in sum_fish[si]:
                   print elem
       print ' ---------------------------------'
       print
       print
       print 'Creating Output root directory if nonexisting'
       
       outdir = os.path.dirname(outroot)
       # create directory if it does not exist
       if outdir=='':
           print('Output root is on working folder')
       elif not os.path.exists(outdir):
           os.makedirs(outdir)
       else: 
           print(str(outdir)+'  exists already')

    if not files[0]:
       if not sum_fish[0][0]:
          print 'NO MATRICES TO WORK WITH!'
          exit()
       else:
          files = []
          #for si in range(summax):
          #    if sum_fish[si][0]:
          #        print("entering extend")
          #        files.extend(sum_fish[si])
          print 'No single fishers to plot, using only the combined ones'
   
    print(fishnames)
    print("files")
    print(files)
    #MOD: too much putput here!
    summing = range(summax)
    
    if derived is not False:
        fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=True)
        for si in range(summax):
            if sum_fish[si][0]:
                print 'NOT HERE'
                summing[si] = fpa.CosmicFish_FisherAnalysis(fisher_path=sum_fish[si], with_derived=True)
    else:
        fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=False)
        for si in range(summax):
            if sum_fish[si][0]:
                print("enter sum pair")
                summing[si] = fpa.CosmicFish_FisherAnalysis(fisher_path=sum_fish[si], with_derived=False)
   


    fishers_temp = fpa.CosmicFish_FisherAnalysis()
    fisher_list = fishers.get_fisher_matrix()
    
    print("Number of fishers before sum: ")
    print(len(fisher_list))
    
    for si in range(summax):
        if sum_fish[si][0]:
            summing_list = summing[si].get_fisher_matrix()
            for fish in summing_list[1:]:
                summing_list[0] = summing_list[0]+fish
            fisher_list.append(summing_list[0])

    print("Number of fishers: ")
    print(len(fisher_list))

    for i in range(len(fisher_list)):
        fisher_list[i].name = fishnames[i]


    fishers_temp.add_fisher_matrix( fisher_list[:] )
    fishers = fishers_temp
 

    #Producing bounds files:
    # get the parameters:
    numbounds     = [ i for i in Config.items( "bounds" ) if "params" in i[0] ]
    numfoms     = [ i for i in Config.items( "FoMs" ) if "params" in i[0] ]
    use_latex     = Config.getboolean( "bounds",'use_latex')
    latex_num_col = Config.getint( "bounds",'latex_num_col')
    
    if len(numfoms)>0:
        if outroot is not None:
            out_file_fom = open(outroot+'_FoMs.txt',"w")
        for key, params in numbounds:
            params = Config.get("FoMs", key).split(",")
            fishers_temp = fishers
            for num, fish in enumerate(fishers_temp.get_fisher_list()):
                out_file_fom.write( 'FoM for ('+', '.join(params)+')\n' )
		print(fish.name)
                out_file_fom.write( str(fish.name)+' \n' )
		fishpars=fish.get_param_names()
                print("all parameter names")
	        print(str(fishpars))
		fishfom=fo.marginalise(fish, params)
                print("FoM for params: ")
		print(str(params))
		fom=fishfom.determinant() 
		fom=np.sqrt(fom)
		print(str(fom))
                out_file_fom.write( str(fom)+' \n' )

        out_file_fom.close()

    if len(numbounds)>0:
    
        if not args.quiet:
            print
            print 'Producing bounds:'
    
        # open the file if wanted:
        if outroot is not None:
            out_file = open(outroot+'_bounds.txt',"w")
    
        
        # do some first printing for Latex:
        if use_latex:
            out_file.write( '\\begin{tabular}{ |'+''.join(['l|' for i in range(latex_num_col) ])+' }\n' )
                
        for key, params in numbounds:
            params = Config.get("bounds", key).split(",")
            
            fishers_temp = fishers
            if eliminate is not False:
               fishers_temp = fishers_temp.reshuffle( params=params )
            elif params is not None:
               fishers_temp = fishers_temp.marginalise( params=params )
            
            for num, fish in enumerate(fishers_temp.get_fisher_list()):
                #output fisher matrices again
                fish.save_to_file(outroot+'__Fisher_Out_'+str(num)+"-"+fish.name, simple_header=True)
		print(fish.name)
		fishpars=fish.get_param_names()
                print("parameter names")
	        print(str(fishpars))
                
                # get the bounds:
                Bounds_68 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.680 ), mode=1 ) )
                Bounds_95 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.950 ), mode=1 ) )
                Bounds_99 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.997 ), mode=1 ) )
                # get the fiducial:
                fiducial = []
                for num, par in enumerate(fish.get_param_fiducial()):
                    fiducial.append( fu.significant_digits( (par, Bounds_68[num]), mode=1 ) )
                # get the names:
                if use_latex:
                    parameter_names_latex = copy.deepcopy(fish.get_param_names_latex())
                else:
                    parameter_names       = copy.deepcopy(fish.get_param_names())
                # do the printing:
                if use_latex:
                    print_table = []
                    for par,fid,bound in zip( parameter_names_latex,fiducial,Bounds_68 ):
                        cent=100
                        psign='\%'
                        if isclose(fid, 0.0, rel_tol=1e-6, abs_tol=1e-6): 
                            fid = 1
                            cent = 1
                            psign=''
                        perc=abs((bound/fid)*cent)
                        print_table.append( '$\sigma_{'+str(par)+'} = '+str('{:4.2f}'.format(perc))+psign+' $' )
                    if len(print_table)%latex_num_col == 0:
                        table_length = len(print_table)/latex_num_col
                    else:
                        table_length = len(print_table)/latex_num_col +1
                    print_table = fu.grouper( table_length, print_table,fillvalue='' )
                    print_table = [ list(i) for i in print_table]      
                    print_table = map(list, zip(*print_table))
                    col_width = [max(len(str(x)) for x in col) for col in zip(*print_table)]
                    
                    if outroot is not None:
                        out_file.write( '\hline'+'\n' )
                        out_file.write( '\multicolumn{'+str(latex_num_col)+'}{|c|}{'+fish.name.replace('_',' ')+'} \\\[1mm]'+'\n' )
                        out_file.write( '\hline'+'\n' )
                        for line in print_table:
                            out_file.write( " " + " & ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + " \\\[1mm]\n" )
                    else:
                        print '\hline'
                        print '\multicolumn{'+str(latex_num_col)+'}{|c|}{'+fish.name.replace('_',' ')+'} \\\[1mm]'
                        print '\hline'
                        for line in print_table:
                            print " " + " & ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + " \\\[1mm]"
                        
                else:
                    # put on top the labels of the columns:
                    parameter_names.insert(0,' Parameter ')
                    fiducial.insert(0, ' fiducial')
                    Bounds_68.insert(0,' 68% C.L.')
                    Bounds_95.insert(0,' 95% C.L.')
                    Bounds_99.insert(0,' 99.7% C.L.')
                    # put a white space:
                    parameter_names.insert(1,' ')
                    fiducial.insert(1, ' ')
                    Bounds_68.insert(1,' ')
                    Bounds_95.insert(1,' ')
                    Bounds_99.insert(1,' ')
                    #
                    print_table = [parameter_names,fiducial,Bounds_68, Bounds_95, Bounds_99]
                    
                    out_file.write( ''.join([ '*' for i in xrange(len('Parameter bounds for the Fisher matrix: '+fish.name)+1)])+'\n' )
                    out_file.write( 'Parameter bounds for the Fisher matrix: '+fish.name+'\n' )
                    out_file.write( ''.join([ '*' for i in xrange(len('Parameter bounds for the Fisher matrix: '+fish.name)+1)])+'\n' )
                    out_file.write( '\n' )
                    print_table = map(list, zip(*print_table))
                    col_width = [max(len(str(x)) for x in col) for col in zip(*print_table)]
                    # print it to file:
                    for line in print_table:
                        out_file.write( "| " + " | ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + " |"+'\n' )
                    out_file.write( '\n' )
                
            if not args.quiet:
               print ' Bounds computed for parameters '+str( fishers_temp.get_parameter_list() )
        
        # finalize the latex part:
        if use_latex:
            out_file.write( '\hline\n' )
            out_file.write( '\end{tabular}' )
        # close the file:
            out_file.close()  
        if not args.quiet:
            print ' Saved results in: ', outroot+'_bounds.txt'
            print 'bounds done!'

    # finalize:
    if not args.quiet:
       print
       print 'It seems everything is done...'
       print 'It was nice working with you. See you soon!'
    # everything's fine, exit without error:
    exit(0)
