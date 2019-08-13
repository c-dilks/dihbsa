hadE = [" "]*2
hadTheta = [" "]*2
hadPhi = [" "]*2
hadFstar = [" "]*2
hadPt = [" "]*2
hadZ = [" "]*2
hadMissMass = [" "]*2
hadXF = [" "]*2
hadEtaCM = [" "]*2
hadEtaBreit = [" "]*2

outStrHad = [" "]*2

#
#
# SEE harut_format_README.txt for columns information
#
#

# three lines for each event: 
# first is pi+ (use s=kP), second is pi- (use s=kM), third is for dihadron
s = 0
kP = 0
kM = 1

with open("xfiles/dihad1000.dis.0000.nrad.dat.evio.hipo2.txt") as iFile, open("xtreeHarut.dat", 'w') as oFile:
    for lineIn in iFile:
        line = lineIn.strip()
        lineCols = line.split()

        #lines with 11 columns are for pi+ and pi-
        if len(lineCols) == 11:
            # s is kP the first time through, kM the second time through
            evnum = lineCols[0]
            hadE[s] = lineCols[1]
            hadTheta[s] = lineCols[2]
            hadPhi[s] = lineCols[3]
            hadFstar[s] = lineCols[4]
            hadPt[s] = lineCols[5]
            hadZ[s] = lineCols[6]
            hadMissMass[s] = lineCols[7]
            hadXF[s] = lineCols[8]
            hadEtaCM[s] = lineCols[9]
            hadEtaBreit[s] = lineCols[10]
            s += 1
            print(line)
        
        #lines with 5 columns are for dihadron kinematics etc.
        if len(lineCols) == 7:
            s = 0 # reset s

            evnum = lineCols[0]
            PhPhiStar = lineCols[1] # actually this is PhiH
            PhPt = lineCols[2]
            phiRT1 = lineCols[3]
            phiRT2 = lineCols[4] # (preferred phiR def)
            # (not sure what lineCols[5] and lineCols[6] are...

            print(line)

            # print output; this is done by appending columns to outStr (with a printout
            # after each append, to easily see what's appended)

            outStr = evnum
            print(outStr)

            for h in range(2):
                outStr = " ".join([outStr,hadE[h],hadTheta[h],hadPhi[h],hadPt[h]])
                outStr = " ".join([outStr,hadXF[h]])
                print(outStr)

            outStr = " ".join([outStr,PhPt,PhPhiStar,phiRT1,phiRT2]);
            print(outStr)

            print "---"
            outStr += "\n"
            oFile.write(outStr)

