#!/usr/bin/python

import sys,math

def massToName(mass):
   elemMass = [1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.180,22.990,24.305,26.982,28.086,30.974,32.066,35.453,39.948,39.098,40.078]
   elemName = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca"]
   for i in range(len(elemMass)):
      if elemMass[i] == mass:
         name = elemName[i]
         break
   return name

optInd = 1
optCount = 0

while optInd < len(sys.argv):
   if sys.argv[optInd] == '-P':
      path = sys.argv[optInd+1]
      optInd = optInd + 2
      optCount = optCount + 1
   elif sys.argv[optInd] == '-n':
      numChunk = eval(sys.argv[optInd+1])
      optInd = optInd + 2
      optCount = optCount + 1
   elif sys.argv[optInd] == '-M':
      numNodes = eval(sys.argv[optInd+1])
      optInd = optInd + 2
      optCount = optCount + 1
   elif sys.argv[optInd] == '-N':
      numPart = eval(sys.argv[optInd+1])
      optInd = optInd + 2
      optCount = optCount + 1
   else:
      if optCount == 0:
         print "Invalid Options!"
         sys.exit(0)
      else:
         break

pathWrite = path + 'traj' + str(numChunk*numNodes) + 'steps.xyz'
auToAng = 0.592
writeXYZ = open(pathWrite,"w")

for i in range(1,numChunk+1):
   fullPath = path + 'chunk' + str(i) + '.traj' 
   readTrajt = open(fullPath,"r")
   coord = []
   mass = []
   while True:
      strRead = readTrajt.readline()
      strArr = strRead.split()
      if len(strArr) == 0:
         break
      mass.append(eval(strArr[2]))
      coord.append(eval(strArr[3]))
      coord.append(eval(strArr[4]))
      coord.append(eval(strArr[5]))

   # WRITE CONFIGURATIONS TO XYZ FILES
   for k in range(numNodes):
      print >> writeXYZ, str(numPart)
      print >> writeXYZ
      for j in range(numPart):
         name = massToName(mass[k*numPart+j])
         print >> writeXYZ, name + ' ' + str(coord[3*k*numPart+3*j]*auToAng) + ' ' + str(coord[3*k*numPart+3*j+1]*auToAng) + ' ' + str(coord[3*k*numPart+3*j+2]*auToAng)

   del coord
   del mass

   readTrajt.close()

writeXYZ.close()
