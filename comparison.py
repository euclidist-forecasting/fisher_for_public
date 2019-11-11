import sys, platform, os
from argparse import ArgumentParser
from time import time
import datetime

#IMPORTING MODULES
sys.path.append('./PlottingScripts')
#sys.path.append('./PlottingScripts/cosmicfishi-pyplots/cosmicfish_pylib')

from appXC import XCcomparison
from appWL import WLcomparison
from appGC import GCcomparisonPostProjection

print ''
print '**************************************************************'
print ' This is the IST Fisher matrix comparison tool.'
print '**************************************************************'
print ''



#GETTING OPTIONS FROM COMMAND LINE
parser = ArgumentParser()
parser.add_argument("-o", "--observables", nargs='+', dest="observables", default=['GC', 'WL', 'XC'],
                  help="Decide on which observables to perform the comparison.")
parser.add_argument("-f","--fisher", nargs='+', type=str, dest="userfisher", default=['TESTmatrix'],
                  help="This is the name of the user provided Fisher matrix. The files must be in the same format of those provided in this repository and use the same name convention (see README).")
args            = parser.parse_args()
userfisher      = args.userfisher
observables     = args.observables

if userfisher == ['TESTmatrix']: 
   print 'The script compares multiple results with respect to their median.'
   print 'Since no external matrix is provided, the script will compare EuclidIST matrices with a TESTmatrix, in order not to have a zero median.'
   print 'This is done only in the XC cases (see Figure 5 of IST paper).'
   print 'To use the other observables, please provide an external Fisher matrix'


empty=0  ## counts if nothing has been done and warns user
if 'GC' in observables: 
    GCcomparisonPostProjection(userfisher)
    empty+=1
if 'WL' in observables: 
    WLcomparison(userfisher)
    empty+=1
if 'XC' in observables: 
    XCcomparison(userfisher)
    empty+=1

if empty==0:
    print("No observables passed with the option -o, no comparison plots were produced")



