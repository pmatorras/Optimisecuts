import os,optparse
from ROOT import TFile
wloc       = os.environ['WWW']

#Optional parameters                                                                                       
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('--year', dest='year', help='which year to open', default='run2')
parser.add_option('--var'   , dest='var'   , help='var to optimise', default='MT2')
parser.add_option('--AN'    , dest='AN'    , help='AN binning', default=False, action='store_true')
parser.add_option('--MT2'   , dest='MT2'   , help='MT2 binning', default=False, action='store_true')
parser.add_option('--Ptm'   , dest='Ptm'   , help='Ptm binning', default=False, action='store_true')
parser.add_option('--all'   , dest='all'   , help='run only few plots', default=False, action='store_true')
parser.add_option('--relerr', dest='relerr', help='do relative error', default=False, action='store_true')

parser.add_option('--sig' , dest='sig' , help='input signal', default='Stop')

(opt, args) = parser.parse_args()



if 'all' in opt.year:
    year = 'run2'
elif opt.year is '0':
    year = '2016'
elif opt.year is '1':
    year = '2017'
elif opt.year is '2':
    year = '2018'
else:
    year = opt.year

varOptim=opt.var
if opt.var not in ['Ptm','MT2']:
    print "wrong variable", opt.var
    exit()

inpfnm= '../rootfiles/plots_HighPtMissOptimisationRegion_'+year+'_SM-'
sig_nm=''

if "Stop" in opt.sig or 'T2tt' in opt.sig:
    inpfnm += 'T2tt_mS-400to700.root'
    sig_nm  = 'T2tt'
    splitnm = 'mS'
    m1min = 400                              
    m1max = 700 
if "SlepSnu" in opt.sig:
    inpfnm += 'TChipmSlepSnu_mC-1to1250.root'
    sig_nm  = 'TChipmSlepSnu'
    splitnm = 'mC'
    m1min = 1                              
    m1max = 1250 

print "INPUT FILE:", inpfnm
Test=not(opt.all)
inpfile = TFile (inpfnm, "READ")
optim   = wloc+"/susy/optimisation/"
folder  = "../Histograms/significance/"+sig_nm+'/'+year
hminMT2 =    0
hmaxMT2 = 1000
hminPtm =    0
hmaxPtm = 2000
nMT2 = 100
nPtm = 100
wPtm = (hmaxPtm-hminPtm)/nPtm
wMT2 = (hmaxMT2-hminMT2)/nMT2
