import os
import sys


def XCcomparison(userfisher):

    print ''
    print 'Doing Cross Correlation comparison'
    print ''


    #if set to True, only show plot, do not save fig
    showonly = False

    names=['EuclidISTF']
    if userfisher != '':
       for newfish in userfisher:
           names.append(newfish)
    colornames=["black","yellow","red","magenta","green","cyan"] 

    colors=' -r '+(' '.join(colornames))
    aliases=" -a "+(" ".join(names))


    casewl = 6

    ######                           ####
    ######  looping over all fishers ####

    for case in range(4):

       #GETTING RID OF THE ABSURD START OF COUNTING FROM ZERO!
       i = case+1

       #Setting the ells
       if i <= 2 :
          suffix='pessimistic'
          dirtolook= 'All_Results/pessimistic/flat'
       else:
          suffix='optimistic' 
          dirtolook= 'All_Results/optimistic/flat'

       if (i == 2) or (i == 4):
          prefix = 'GCsp_'
       else:
          prefix = ''

       filebase = prefix+'GCph_WL_XC_w0wa_flat_'+suffix

       if i == 2:
          parlist = ['Omegam','Omegab','w0','wa','h','ns','sigma8','aIA','etaIA','betaIA','b1','b2','b3','b4','b5']
       else:
          parlist = ['Omegam','Omegab','w0','wa','h','ns','sigma8','aIA','etaIA','betaIA','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10']


       if os.path.isfile(dirtolook+'/'+newfish+'_'+filebase+'.txt'):
          print('using',parlist)

          notefile = filebase

          print('working on case: '+filebase)

          fishernames=[namm+'_'+filebase+'.txt' for namm in names]

          filenames=[os.path.join(dirtolook, onefile) for onefile in fishernames]
        
          fishfiles=" ".join(filenames)
       
          dirtosave = dirtolook+'/results_XC/'

          plerr=" -e -m "

          plcase=" -c "+str(casewl)

          noteplot=" -n "+notefile

          plparams=" -p "+(" ".join(parlist))

          runscript="python PlottingScripts/plfmat.py --notstamp_outdir -o "+dirtosave+aliases+colors+plcase+noteplot+plparams+plerr+fishfiles

          os.system(runscript)

          print ''
          print 'Cross Correlation comparison case {} done'.format(filebase)
          print ''       
       else:
          print 'No Fisher matrix to compare with'
          print 'Comparison case {} skipped'.format(filebase)

    return None
