""" genhead.py
	Generate a header if somebody forgot.
	18/03/2016 Markovic
"""

import fishlib as l
import argparse, sys
import fileinput # http://stackoverflow.com/questions/3162314/add-headers-to-a-file

def resave(filename, headerline):
	with open(filename, 'r') as f: filestring = f.read()
	with open(filename, 'w') as f: f.write(headerline +'\n' + filestring)

def genhead(pzs,zs,zorder,exclude,ps,shape_first=True):
	if shape_first:
		header = '# ' + ' '.join([ ' '.join(l.mixpars(pzs,zs,zorder,exclude)), ' '.join(ps) ])
	else:
		header = '# ' + ' '.join([ ' '.join(ps), ' '.join(l.mixpars(pzs,zs,zorder,exclude)) ])
	return header

if __name__=='__main__':

	parser = argparse.ArgumentParser(description="Generate a header to be compied into a file")
	parser.add_argument("-p", "--ps", default=l.PARS, nargs='+', help='list of redshift-independent parameters')
	parser.add_argument("-d", "--pzs", default=l.ZPARS, choices=l.ZPARS, nargs='+', help='list of redshift-dependent parameters')
	parser.add_argument("-n", "--nlin", action = 'store_true', default=False, help='includes nonlinear parameters sigp and sigv if flag is present')
	parser.add_argument("-l", "--dolens", action='store_true', default=False, help='do the full set of lensing parameters instead')
	parser.add_argument('-z', '--zs', type=float, nargs='+', default=l.ZS, help='list of redshift')
	parser.add_argument("-o", "--zorder", action='store_false', default=True, help='order: default is par1_z1, par2_z1...')
	parser.add_argument('-m', '--fidmod', type=str, nargs='+', help='list of fiducial values (must match the order of the parameters: can test without this first)')
	parser.add_argument("-e", "--exclude", nargs='*', default=[], help='optionally exclude one or more of ' + ' '.join(l.PARS))
	parser.add_argument("-f", "--file", nargs='*', default=[], help='add header line to files directly')
	args = parser.parse_args()

	# Get the order where the shape parameters appear at the command line
	if '-p' in sys.argv:
		pind = sys.argv.index('-p')
	elif '--ps' in sys.argv:
		pind = sys.argv.index('--ps')
	else:
		pind = 0

	# Get the order where the z-dependent parameters appear at the command line
	if '-d' in sys.argv:
		dind = sys.argv.index('-d')
	elif '--pzs' in sys.argv:
		dind = sys.argv.index('--pzs')
	else:
		dind = 0
                args.pzs = []

	# Include non-linear parameters into the shape list as they behave the same
	if args.nlin: args.ps = args.ps+l.NLPARS

	# Add lensing nuisance if requested and set the z-dependent set to empty
	if args.dolens: 
		args.pzs = []
		args.ps += l.WLNUISANCE

	# Create the header string
	header = genhead(args.pzs,args.zs,args.zorder,args.exclude,args.ps,pind>dind)

	# Add fiducial values if they are given
	if args.fidmod is not None:
		header += '\n# ' + ' '.join(args.fidmod)

	# Optionally save into the file header if filename given
	if len(args.file)>0:
		for file in args.file: 
			resave(file, header)
			print "Header: \n" + header + "\nAdded to file:  " + file + "."
	else:
		print header
