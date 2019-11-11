import os
import sys



def WLcomparison(userfisher):

    print ''
    print 'Doing Weak Lensing comparison'
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
    ######  looping over all cases   ####

    for case in range(4):

       #GETTING RID OF THE ABSURD START OF COUNTING FROM ZERO!
       i = case+1

       #Setting the ells
       if i <= 2 :
          suffix='pessimistic'
       else:
          suffix='optimistic'

       if (i == 2) or (i == 4):
          dirtolook= 'All_Results/'+suffix+'/non-flat'

          model = 'w0wa_nonflat'
          parlist = ['Omegam','Omegade','Omegab','w0','wa','h','ns','sigma8','aIA','etaIA','betaIA']
          casewl = 7
       else:
          dirtolook= 'All_Results/'+suffix+'/flat'

          model = 'w0wa_flat'
          parlist = ['Omegam','Omegab','w0','wa','h','ns','sigma8','aIA','etaIA','betaIA']
          casewl = 6

       filebase = 'WL_'+model+'_'+suffix



       if os.path.isfile(dirtolook+'/'+newfish+'_'+filebase+'.txt'):
          print('using',parlist)

          notefile = filebase

          print('working on case: '+filebase)

          fishernames=[namm+'_'+filebase+'.txt' for namm in names]

          filenames=[os.path.join(dirtolook, onefile) for onefile in fishernames]

          fishfiles=" ".join(filenames)

          dirtosave = dirtolook+'/results_WL/'

          plerr=" -e -m "

          plcase=" -c "+str(casewl)

          noteplot=" -n "+notefile

          plparams=" -p "+(" ".join(parlist))

          runscript="python PlottingScripts/plfmat.py --notstamp_outdir -o "+dirtosave+aliases+colors+plcase+noteplot+plparams+plerr+fishfiles

          os.system(runscript)

          print ''
          print 'WL comparison case {} done'.format(filebase)
          print ''
       else:
          print 'No Fisher matrix to compare with'
          print 'WL comparison case {} skipped'.format(filebase)



    return None
