""" fishlib.py
    Utilities for plotting the results of the code comparison:
    Markovic, 18/03/2016
    Major updates:
    Casas, 20/07/2017
"""

import matplotlib.pyplot as plt
import argparse
import numpy as np
import os.path as op
from scipy.sparse.linalg import inv
import time, datetime
import subprocess as sp
import sys

# Precision for Fisher matrix
ROUNDTO = 8

# Testing cases
CASES = [1,
         2,
         3,
         4, # All redshift-dependent parameters allowed to vary
         5, # Plot 4 shape parameters and 5 redshift-dependent parameters
                 6, # Final matrices, flat case, LCDM, w0, wa
                 7, # Final matrices, non-flat case, LCDM, w0, wa
                 8, # Final matrices, flat case, LCDM+gamma, w0, wa
                 9, # Final matrices, non-flat case, LCDM+gamma, w0, wa
         10, # 4 shape, 2 nonlinear
         11, # 4 shape, 2 nonlinear and 5xN z-dependent parameters
         12] #XC case with redshift cut

# Parameters (order assumed)
# Betta ZPARS = ['lnfs8', 'lnbs8', 'lnDa', 'lnH', 'lnPs']
PARS = ['wm', 'wb' ,'h', 'ns'] # default backwards-compatible set of parameters
shapePARS = ['wm', 'wb' ,'h', 'ns'] # shape parameters (redshift-independent)
flatPARS = ['Omegam','Omegab','w0','wa','h','ns','sigma8']
nonflatPARS = ['Omegab','h','Omegam','ns','Omegad','w0','wa','sigma8']
flatgammaPARS = ['Omegam','Omegab','w0','wa','h','ns','sigma8', 'gamma']
nonflatgammaPARS = ['Omegab','h','Omegam','ns','Omegad','w0','wa','sigma8', 'gamma']
ZPARS = ['lnDa', 'lnH', 'lnfs8', 'lnbs8', 'Ps'] # redshift-dependent parameters
#ZS = np.arange(0.95, 1.85, 0.1)
ZS = np.array([1.0,1.2,1.4,1.65])
NLPARS = ['sp','sv']
WLNUISANCE = ["aIA", "etaIA", "betaIA"]
biaspars = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10']
cutbias = ['b1', 'b2', 'b3', 'b4', 'b5']

FIDUS_BASELINE = {
      '10e9As': 2.12605,
      'omegab': 0.022445,
      'Omegab': 0.05,
      'Omegad': 0.68,
      'w0': -1.0,
      'wa': 0.0,
      'ns': 0.96,
      'h': 0.67,
      'omegam': 0.143648,
      'wm': 0.143648,
      'wb': 0.022445,
      'Omegam': 0.32,
      'sigma8': 0.815583,
      'gamma': 0.55,
      'aIA': 1.72,
      'etaIA': -0.41,
      'betaIA': 2.17
       }


TeXLabs = {
      '10e9As': r'$10^9A_s$',
      'omegab': r'$\omega_{\rm{b},\!0}$',
      'Omegab': r'$\Omega_{\rm{b},\!0}$',
      'omegac': r'$\omega_{\rm{c},\!0}$',
      'Omegac': r'$\Omega_{\rm{c},\!0}$',
      'Omegad': r'$\Omega_{\rm{DE},\!0}$',
      'omegam': r'$\omega_{\rm{m},\!0}$',
      'Omegam': r'$\Omega_{\rm{m},\!0}$',
      'wm': r'$\omega_{\rm{m},\!0}$',
      'wb': r'$\omega_{\rm{b},\!0}$',
      'ns': r'$n_{\rm{s}}$',
      'w0': r'$w_0$',
      'wa': r'$w_a$',
      'h': r'$h$',
      'sigma8': r'$\sigma_{8}$',
      'gamma': r'$\gamma$',
      'sp':r'$\sigma_{p}$',
      'sv':r'$\sigma_{v}$',
      'aIA':r'$a_{\rm IA}$',
      'etaIA':r'$\eta_{\rm IA}$',
      'betaIA':r'$\beta_{\rm IA}$',
      'lnDa': '\ln D_{A}',
      'lnH': '\ln H',
      'lnfs8': '\ln f\sigma_{8}',
      'lnbs8': '\ln b\sigma_{8}',
      'Ps': 'P_{s}'
       }
# Just a list of colors
COLOURS = 'bgrcmykbgrcmykbgrcmykbgrcmyk'

# Path to this file
__PATH__ = op.dirname(__file__)


def reverse_timestamp(tstampstr):
    t = (2016,1,1,0,0,0,0,0,0)
    sec = int(tstampstr, 16)
    secsince = sec + int(time.mktime(t))
    timestr = time.ctime(secsince)
    print("Time stamp was set at: "+timestr)
        

def get_timestamp():
    t = (2016,1,1,0,0,0,0,0,0)
    return str(hex(int(time.time() - time.mktime(t))))[2:]

def reorder(incols, inlist, outcols):
    if len(incols) != len(inlist): raise Exception("The number of in columns you specified ("+str(len(incols))+") dones't match the length of the list you gave me ("+str(len(inlist))+")!")
    if len(incols) < len(outcols): raise Exception("You are asking me to pick " + str(len(outcols)) + " parameters out of " + str(len(incols)) + "!")
    order = []
    for col in outcols:
        order.append(incols.index(col))
    return list[np.array(inlist)[order]]

def mreorder(incols, inarr, outcols):
    if len(incols) != inarr.shape[1]: raise Exception("The number of in columns you specified ("+str(len(incols))+") doesn't match the number of columns in the array you gave me ("+str(inarr.shape[1])+")!")
    if len(incols) < len(outcols): raise Exception("You are asking me to pick " + str(len(outcols)) + " parameters out of " + str(len(incols)) + "!")
    order = []
    for col in outcols:
        order.append(incols.index(col))
    return inarr[:,order][order,:] 

def marginalise(incols, data, keep=None):
    data = mreorder(incols, np.linalg.inv(data), keep)
    return np.linalg.inv(data)

def mixpars(ps=ZPARS,zs=ZS,zorder=True,exclude=[]):
    if zorder:
        return parszd(ps,zs,exclude)
    else:
        return zdpars(ps,zs,exclude)

def zdpars(ps=ZPARS,zs=ZS,exclude=[]):
    parlist = []
    for p in ps:
        if p not in exclude and p in PARS:
            parlist.append(p)
            continue
        for z in zs:
            if p not in exclude and p in ZPARS: 
                parlist.append(p+'_'+str('{0:.2f}'.format(z)))
    return parlist

def parszd(ps=ZPARS,zs=ZS,exclude=[]):
    parlist = []
    for z in zs:
        for p in ps:
            if p not in exclude and p in ZPARS: 
                parlist.append(p+'_'+str('{0:.2f}'.format(z)))
            elif p in PARS and p not in parlist:
                parlist.append(p)
    return parlist

def whichpz(pars):
    p=[]
    z=[]
    for par in pars:
        par = par.split('_')
        if par[0] not in p: p.append(par[0])
        zt = float(par[1])
        if zt not in z: z.append(zt)

    # Do a check
    missing = []
    extra = []
    checkpars = zdpars(p,z)
    for par in checkpars:
        if par not in pars: missing.append(par)
    for par in pars:
        if par not in checkpars: extra.append(par)
    if len(missing)!=0 or len(extra)!=0:
        raise Exception("Your file is missing the parameter: " + str(missing) + "\nand it contains the following unknown parameters: " + str(extra) + "!")

    return p, z

def selectpars(inpars=ZPARS, keep=['lnDa', 'lnH'], dropped=False):
    outpars=[]
    drop=[]
    for par in inpars:
        appended = False
        for k in keep:
            if k in par: 
                outpars.append(par)
                appended = True
        if not appended:
            drop.append(par)
    if dropped:
        return outpars, drop
    else:
        return outpars

def shorten(instring,length=20):
    """ Shorten string by adding ... for long legend entries. """
    if len(instring)>length:
        instring = instring[:length//2] + '...' + instring[-length//2:]
    return instring

def readDomenico(filestr):
    " Domenico''s format."
    rows = filestr.split('{')
    count=0
    data = []
    for row in rows:
        line = []
        for el in row.replace('{','').replace('}','').split(','):
            el = el.strip()
            if el == '': continue
            line.append(float(el.replace('\n','').strip().replace('*^','E')))
        if len(line)>0: data.append(line)
    return data

def correct_lnPs(parlist, filename=False):
    # Correct any lnPs to Ps, since lnPs makes no sense
    Pchanged = False
    newpars = []
    for par in parlist:
        if 'lnPs' in par: 
            par = par.replace('lnPs', 'Ps')
            Pchanged = True
        newpars.append(par)
    
    if Pchanged:
        printstring = "I'm assuming you meant Ps, not lnPs"
        if filename: 
            printstring += " in "+filename+"!"
        else:
            printstring += "!"
        print printstring

    return newpars

def readhead(filename):
    header = None
    with open(filename) as f:
        filestr = f.read().strip()

    if '#' in filestr:
        header = filestr.split('\n')[0]
        parlist = header.replace('#','').strip().strip('"').split()
    else:
        raise Exception("No header in " + filename + "!")

    print parlist
    parlist = correct_lnPs(parlist, filename)
    print parlist

    return parlist

def readfish(filename):
    header = None
    with open(filename) as f:
        filestr = f.read().strip()

    if '#' in filestr:
        header = filestr.split('\n')[0]
        filestr = filestr.replace(header,'')
        parlist = header.replace('#','').strip().strip('"').split()
    else:
        raise Exception("No header in " + filename + "!")

    # Correct any lnPs to Ps, since lnPs makes no sense
    parlist = correct_lnPs(parlist, filename=False)

    if '{{' in filestr:
        data = readDomenico(filestr)
    else:
        data = np.genfromtxt(filename)

    return np.matrix(data), parlist

def readerrs(filename, parlist=zdpars()):

    with open(filename) as f:
        headstr = f.readline().strip()

    if '#' in headstr:
        header = filestr.split('\n')[0]
        filestr = filestr.replace(header,'')

    data = np.genfromtxt(filename, names=True)

    return data

def readpks(filename, parlist=zdpars()):

    data = np.genfromtxt(filename, names=True)

    return data


def recognizecase(casenum):
    if casenum == 5:
        paras = shapePARS+ZPARS
    if casenum == 6:
        paras = flatPARS+WLNUISANCE+biaspars
    if casenum == 7:
        paras = nonflatPARS
    if casenum == 8:
        paras = flatgammaPARS
    if casenum == 9:
        paras = nonflatgammaPARS
    if casenum == 10:
        paras = shapePARS+NLPARS
    if casenum == 11:
        paras = shapePARS+NLPARS+ZPARS
    if casenum == 12:
        paras = flatPARS+WLNUISANCE+cutbias
    if casenum == 60:
        paras = flatPARS+WLNUISANCE
    if casenum == 70:
        paras = nonflatPARS+WLNUISANCE
    if casenum == 80:
        paras = flatgammaPARS+WLNUISANCE
    if casenum == 90:
        paras = nonflatgammaPARS+WLNUISANCE
    return paras
    


    # A class that contains the Git environment at the time of it's initialisation.
    # Currently it uses the subprocess module to speak to Git through the system.
    #   Ideally some day it would use the GitPython module.
class GitEnv(object):

    # Optional input, /path/to/.git
    def __init__(self, git_dir=None):
        self.git_dir = git_dir
        self.hash, self.author, self.date = [str(s) for s in self.get_commit()]
        self.url = str(self.get_remote())
        self.branch = str(self.get_branch())
        self.repo = str(self.get_repo())
        self.printstart = ''
    # Also, should have an if that gives out the name of the parent folder + the
    # date and time in the case that it is NOT A GIT REPO!
    
    def __str__(self):
        startline = self.printstart
        as_string = startline + "This was generated by code from the Git repository at:"
        as_string += "\n" + startline + "\t " + self.url + ","
        as_string += "\n" + startline + "\t on the " + self.branch + " branch,"
        as_string += "\n" + startline + "\t with commit: " + self.hash
        as_string += "\n" + startline + "\t\t from " + self.date + ", "
        as_string += "\n" + startline + "\t\t by " + self.author + "."
        return as_string

    def set_print(self, startline):
        self.printstart = startline
    
    def get_git_cmd(self, args=[]):
        cmd = ['git']
        if self.git_dir != None:
            cmd.append('--git-dir')
            cmd.append(self.git_dir)
        for one in args:
            cmd.append(one)

        return cmd

    def get_hash(self, nochar=7, sep=''):
        return sep+self.hash[0:nochar]+sep
    
    # Get the hash, author and date of the most recent commit of the current repo.
    def get_commit(self):
        cmd = sp.Popen(self.get_git_cmd(['log', '-n','1']), stdout=sp.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        newlist=[]
        for entry in cmd_out.strip().split('\n'):       
            if entry=='': continue
            entry = entry.split(' ')
            # This is a hack, should use a dict so can be sure what we are reading in:
            if 'commit' in entry[0] or 'Author' in entry[0] or 'Date' in entry[0]:
                newlist.append(' '.join(entry[1:]).strip())
        return newlist

    # At the moment this only gives the first url in what git returns.
    # Eventually it'd be nice if you could get the origin url, the fetch...
    def get_remote(self):
        cmd = sp.Popen(self.get_git_cmd(['remote', '-v']), stdout=sp.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        if bool(cmd_out):
            try:
                return cmd_out.strip().split('https://')[1].split(' ')[0]
            except IndexError:
                ssh_url = cmd_out.strip().split('git@')[1].split(' ')[0]
                return ssh_url.replace(':','/')
        else:
            return 'no remote repo'

    def get_branch(self):
        cmd = sp.Popen(self.get_git_cmd(['branch']), stdout=sp.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        branches = cmd_out.strip().splitlines()
        for branch in branches:
            if '*' in branch:
                return branch.replace('*','').strip()

    def get_repo(self):
        cmd = sp.Popen(self.get_git_cmd(['rev-parse','--show-toplevel']), stdout=sp.PIPE)
        cmd_out, cmd_err = cmd.communicate()
        return cmd_out.strip().split('/')[-1]

def savelog(inargs, filename='./log.txt', info=''):
    git = GitEnv()
    git.set_print('')

    info = info + '\n' + str(datetime.datetime.now()) + ' UTC'
    info = info + '\n' + ',\n'.join(str(inargs).replace('Namespace(','').replace(')','').split(', '))
    info = info + '\n' + str(git)

    with open(filename,'w') as f:
        f.write(info)

    print 'Printed log to ' + filename + '.'

    return 0

if __name__=='__main__':
    print whichpz(readhead(sys.argv[1]))
