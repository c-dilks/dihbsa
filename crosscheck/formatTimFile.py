import re
import math

eleP = 0
pipP = 0
pimP = 0

def getVal(col):
  return col.split(": ")[-1]

with open("xfiles/hayward_cross_check_0001.txt") as iFile, open("xtree.dat", 'w') as oFile:
  for lineIn in iFile:
    line = lineIn.strip()
    lineCols = line.split(",")

    if "Event:" in line:
      evNum = getVal(lineCols[0])
      #print(line)
      #print(evNum)
      inStr = lineIn

    if "Q2:" in line:
      Q2 = getVal(lineCols[0])
      W = getVal(lineCols[1])
      x = getVal(lineCols[2])
      y = getVal(lineCols[3]).split()[0]
      #print(line)
      #print(Q2,W,x,y)
      inStr += lineIn

    if "PID:" in line:
      lineCols = line.split(" ")
      pid = int(getVal(lineCols[1]))
      mom = float(re.sub(",","",getVal(lineCols[3])))
      if pid == 11:
        eleP = mom if mom>eleP else eleP
      elif pid == 211:
        pipP = mom if mom>pipP else pipP
      elif pid == -211:
        pimP = mom if mom>pimP else pimP
      inStr += lineIn

    if "pair mass:" in line:
      Mh = getVal(lineCols[0])
      pT = getVal(lineCols[1])
      xF = getVal(lineCols[2])
      theta = getVal(lineCols[3])
      phiR = getVal(lineCols[4])
      phiH = getVal(lineCols[5])
      #print(line)
      #print(Mh,pT,xF,theta,phiR,phiH)
      inStr += lineIn

      eleE = str( math.sqrt( eleP**2 + 0.000511**2 ) )
      pipE = str( math.sqrt( pipP**2 + 0.139571**2 ) )
      pimE = str( math.sqrt( pimP**2 + 0.139571**2 ) )

      outStr = " ".join([evNum,eleE,pipE,pimE,Q2,W,x,y,Mh,pT,xF,theta,phiR,phiH])
      outStr += "\n"
      print(inStr)
      print(outStr)
      print("-----")
      oFile.write(outStr)

      eleP = 0
      pipP = 0
      pimP = 0


        

      

