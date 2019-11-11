import os
import sys
import fileinput
import shutil

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

def replaceInFile(filename,searchExp,replaceExp):
    for line in fileinput.input(filename, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)

thirteen='             '
newline = '\n'
fishpos='sum_fish4:    '   #for configuration ini file

print("running python cosmicfish plotting code")
if 'PlottingScripts/' in os.getcwd():
    basefolder=''
else:
    basefolder='PlottingScripts/'

inifile_all='ellipses-FoMs-optimistic-flat-w0waCDM.ini'
inifile_xc='xc_ellipses-FoMs-optimistic-flat-w0waCDM.ini'
inifile_gc='gc_ellipses-FoMs-optimistic-flat-w0waCDM.ini'
inifile_wl='wl_ellipses-FoMs-optimistic-flat-w0waCDM.ini'

def EllipsesPlot(userfisher='', case='All', name='', color=''):

    if case=='All':
        appy='ellipses_calc.py'
        inifile=inifile_all
    if case=='XC':
        appy='xc_ellipses_calc.py'
        inifile=inifile_xc
    if case=='GC':
        appy='xc_ellipses_calc.py'
        inifile=inifile_gc
    if case=='WL':
        appy='xc_ellipses_calc.py'
        inifile=inifile_wl

    if userfisher != '':
        print("user fisher added")
        print("Relative path and filename to user Fisher matrix: ")
        print(userfisher)
        print("Name of Fisher : ")
        if name=='':
            name = 'USER'
        print(name)
        print("User color : ")
        if color=='':
            color = 'black'
        print(color)
        userinifile = 'USER_'+inifile
        print(" Creating new USER ini file: ")
        shutil.copy2(basefolder+inifile, basefolder+userinifile)
        print(basefolder+userinifile)

        replaceInFile(basefolder+userinifile,"USERPATH",fishpos+userfisher+newline)
        replaceInFile(basefolder+userinifile,"USERFISHERNAME",thirteen+name+newline)
        replaceInFile(basefolder+userinifile,"USERCOLOR",thirteen+color+newline)      
        
        inifile_full=basefolder+userinifile
    else:
        inifile_full = basefolder+inifile
    
    
    runscript="python "+basefolder+'cosmicfish-pyplots/apps/'+appy+' '+inifile_full

    os.system(runscript)
    print("ellipses and fom results done")
    return 1

if __name__ == '__main__':
    EllipsesPlot()
