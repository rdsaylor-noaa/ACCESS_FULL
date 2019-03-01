#!/usr/bin/env python
#
# fsplit.py
#
# Reads ACCESS output file and splits into separate 
# single species files for plotting with ncl
#
#=====================================================================================!
#                                                                                     !
#     Program:      ACCESS                                                            !
#                   Atmospheric Chemistry and Canopy Exchange Simulation System       !
#                                                                                     !
#     Version:      3.1.0                                                             !
#                                                                                     !
#     Last Update:  February 2019                                                     !
#                                                                                     !
#     Contact:      Rick D. Saylor, PhD                                               !
#                   Physical Scientist                                                !
#                   U. S. Department of Commerce                                      !
#                   National Oceanic and Atmospheric Administration                   !
#                   Air Resources Laboratory                                          !
#                   Atmospheric Turbulence and Diffusion Division                     !
#                   456 S. Illinois Ave                                               !
#                   Oak Ridge, TN 37830                                               !
#                   email: Rick.Saylor@noaa.gov                                       !
#                                                                                     !
#=====================================================================================!
#                                                                                     !
import os
import sys

# usage
def usage():
    print("usage: $s SIMNAME TYPE" % os.path.basename(sys.argv[0]))

def main(argv=None):
    if argv is None:
        argv = sys.argv
 
    # enforce proper usage of 2 and only 2 arguments
    if len(argv) != 3:
        usage()
        return 2

    # get filename and attach filehandle
    simname = argv[1]
    type = argv[2]

    if type == 'budget':

        f1 = simname+'/budget/bcn.out'
        print("Reading file: %s" % f1)
        fh1 = open(f1)
        lines1 = fh1.readlines()
        for l1 in lines1:
            tokens = l1.split()
            if l1[0] != ' ':
                fh1.close()
                of1 = simname+'/budget/'+tokens[0]+'_bcn.dat' 
                fh1 = open(of1, 'w')
                continue
            fh1.write(l1)
        fh1.close()

        f2 = simname+'/budget/bch.out'
        print("Reading file: %s" % f2)
        fh2 = open(f2)
        lines2 = fh2.readlines()
        for l2 in lines2:
            tokens = l2.split()
            if l2[0] != ' ':
                fh2.close()
                of2 = simname+'/budget/'+tokens[0]+'_bch.dat' 
                fh2 = open(of2, 'w')
                continue
            fh2.write(l2)
        fh2.close()

        f3 = simname+'/budget/bdp.out'
        print("Reading file: %s" % f3)
        fh3 = open(f3)
        lines3 = fh3.readlines()
        for l3 in lines3:
            tokens = l3.split()
            if l3[0] != ' ':
                fh3.close()
                of3 = simname+'/budget/'+tokens[0]+'_bdp.dat' 
                fh3 = open(of3, 'w')
                continue
            fh3.write(l3)
        fh3.close()

        f4 = simname+'/budget/bem.out'
        print("Reading file: %s" % f4)
        fh4 = open(f4)
        lines4 = fh4.readlines()
        for l4 in lines4:
            tokens = l4.split()
            if l4[0] != ' ':
                fh4.close()
                of4 = simname+'/budget/'+tokens[0]+'_bem.dat' 
                fh4 = open(of4, 'w')
                continue
            fh4.write(l4)
        fh4.close()

        f5 = simname+'/budget/bvt.out'
        print("Reading file: %s" % f5)
        fh5 = open(f5)
        lines5 = fh5.readlines()
        for l5 in lines5:
            tokens = l5.split()
            if l5[0] != ' ':
                fh5.close()
                of5 = simname+'/budget/'+tokens[0]+'_bvt.dat' 
                fh5 = open(of5, 'w')
                continue
            fh5.write(l5)
        fh5.close()

    elif type == 'vs':
        fname = simname+'/r/'+type+'.out'
        print("Reading file: %s" % fname)
        fh = open(fname)

        # read lines from fname
        lines = fh.readlines()
        fh.close()

        # loop over each line, write each species data to separate file
        for line in lines:
            tokens = line.split()
            if (len(tokens) == 1):
                ofname = simname+'/r/'+tokens[0]+'_vs.dat'
                fh = open(ofname, 'w') 
                continue
            if (tokens[0] == 'z(m)'):
                hrs = tokens[1:]
            else:
                vs  = tokens[1:]
                fh.write('    hr       vs(cm/s)\n')
                for ihr in range(len(hrs)):
                   fh.write('%6s       %s\n' % (hrs[ihr], vs[ihr]))
                fh.close() 
                
    elif type == 'rsoill':
        fname = simname+'/r/'+type+'.out'
        print("Reading file: %s" % fname)
        fh = open(fname)

        # read lines from fname
        lines = fh.readlines()
        fh.close()

        # loop over each line, write each species data to separate file
        for line in lines:
            tokens = line.split()
            if (len(tokens) == 1):
                ofname = simname+'/r/'+tokens[0]+'_rsoill.dat'
                fh = open(ofname, 'w') 
                continue
            if (tokens[0] == 'z(m)'):
                hrs = tokens[1:]
            else:
                rsoil  = tokens[1:]
                fh.write('    hr       rsoill(s/cm)\n')
                for ihr in range(len(hrs)):
                   fh.write('%6s       %s\n' % (hrs[ihr], rsoil[ihr]))
                fh.close() 
                
    elif type == 'rcanopy':
        f1 = simname+'/r/rb.out'
        print("Reading file: %s" % f1)
        fh1 = open(f1)
        lines1 = fh1.readlines()
        for l1 in lines1:
            tokens = l1.split()
            if l1[0] != ' ':
                fh1.close()
                of1 = simname+'/r/'+tokens[0]+'_rb.dat' 
                fh1 = open(of1, 'w')
                continue
            fh1.write(l1)
        fh1.close()

        f2 = simname+'/r/rs.out'
        print("Reading file: %s" % f2)
        fh2 = open(f2)
        lines2 = fh2.readlines()
        for l2 in lines2:
            tokens = l2.split()
            if l2[0] != ' ':
                fh2.close()
                of2 = simname+'/r/'+tokens[0]+'_rs.dat' 
                fh2 = open(of2, 'w')
                continue
            fh2.write(l2)
        fh2.close()

        f3 = simname+'/r/rc.out'
        print("Reading file: %s" % f3)
        fh3 = open(f3)
        lines3 = fh3.readlines()
        for l3 in lines3:
            tokens = l3.split()
            if l3[0] != ' ':
                fh3.close()
                of3 = simname+'/r/'+tokens[0]+'_rc.dat' 
                fh3 = open(of3, 'w')
                continue
            fh3.write(l3)
        fh3.close()

        f4 = simname+'/r/rm.out'
        print("Reading file: %s" % f4)
        fh4 = open(f4)
        lines4 = fh4.readlines()
        for l4 in lines4:
            tokens = l4.split()
            if l4[0] != ' ':
                fh4.close()
                of4 = simname+'/r/'+tokens[0]+'_rm.dat' 
                fh4 = open(of4, 'w')
                continue
            fh4.write(l4)
        fh4.close()

    else:
        fname = simname+'/'+type+'/'+type+'.out'
        print("Reading file: %s" % fname)
        fh = open(fname)

        # read lines from fname
        lines = fh.readlines()

        # loop over each line, write each species data to separate file
        for line in lines:
            tokens = line.split()
            if line[0] != ' ': 
                fh.close() 
                ofname = simname+'/'+type+'/'+tokens[0]+'.dat'
                fh = open(ofname, 'w')  
                continue
            fh.write(line)     
        fh.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
