import optparse
import sys, os
cmsenv=' eval `scramv1 runtime -sh` '
optim='/afs/cern.ch/work/p/pmatorra/private/CMSSW_10_2_14/src/Optimisecuts/'

if __name__ == '__main__':
    doDC=False
    if len(sys.argv)>4:
        print 'Please, specify Sample, number of events and file location, in that order'
        sys.exit()
    sample="T2tt"
    nSam  ="1000"
    sfile ="no"
    if(len(sys.argv)>1):
       sample=sys.argv[1]
       if(len(sys.argv)>2): 
           nSam=sys.argv[2]
           if(len(sys.argv)>3): 
               sfile=sys.argv[3]
    
       
    runcommand="python loop.py --sample="+sample+" --nsample="+nSam+" --samplefile="+sfile
    os.system("cd "+optim+" ; "+cmsenv+" ; "+runcommand)
