import sys, platform, os
from argparse import ArgumentParser
from time import time
import datetime

#IMPORTING MODULES
sys.path.append('./PlottingScripts')
#sys.path.append('./PlottingScripts/cosmicfishi-pyplots/cosmicfish_pylib')

from calcEllipsesFoMs import EllipsesPlot

print ''
print '**************************************************************'
print ' This is the IST Fisher matrix ellipses plotting tool.'
print '**************************************************************'
print ''



#GETTING OPTIONS FROM COMMAND LINE
parser = ArgumentParser()
parser.add_argument("-o", "--observables", nargs='+', dest="observables", default=['All'],
                  help="Decide which observables to plot.")
parser.add_argument("-f","--fisher", nargs='?', type=str, dest="userfisher", default='',
                  help="This is the name of the user provided Fisher matrix. The files must be in the same format of those provided in this repository and use the same name convention (see README).")

parser.add_argument("-c","--color", nargs='?', type=str, dest="usercolor", default='',
                  help="This is the color for the elliptical contour of the Fisher matrix provided by the user.")

parser.add_argument("-n","--name", nargs='?', type=str, dest="username", default='',
                  help="This is the name for the Fisher matrix provided by the user.")
args            = parser.parse_args()
userfisher      = args.userfisher
observables     = args.observables
usercolor       = args.usercolor
username        = args.username

empty=0  ## counts if nothing has been done and warns user

if 'All' in observables:
    print("**** Ellipses for all combined probes****")
    EllipsesPlot(userfisher, case='All', color=usercolor, name=username)
    empty+=1

if 'GC' in observables:
    print("**** Ellipses for GC probes****")
    EllipsesPlot(userfisher, case='GC', color=usercolor, name=username)
    empty+=1

if 'WL' in observables: 
    print("**** Ellipses for WL probes****")
    EllipsesPlot(userfisher, case='WL', color=usercolor, name=username)
    empty+=1

if 'XC' in observables: 
    print("**** Ellipses for XC probes****")
    EllipsesPlot(userfisher, case='XC', color=usercolor, name=username)
    empty+=1

if empty==0:
    print("No observables passed with the option -o, no comparison plots were produced")



