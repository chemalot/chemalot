#!/usr/bin/env python
#Alberto
from __future__ import print_function, division
import os
import glob
import argparse
from subprocess import Popen, PIPE
from argparse import RawTextHelpFormatter
import sys
import re
from textwrap import dedent
from subprocess import call


def warn(*objs):
   print(*objs, file=sys.stderr)


class FileLineWrapper(object):
   """File reader with line numbers and readRequiredLine"""

   def __init__(self, f):
      self.f = f
      self.line = 0

   def close(self):
      return self.f.close()

   def readline(self):
      self.line += 1
      return self.f.readline()

   def readRequiredLine(self):
      line = self.readline()
      if not line:
         raise IOError("Unexpected end of file found: %d" % self.line)
      return line

   def __enter__(self):
      self.f.__enter__()
      return self

   def __exit__(self, tp, value, traceback):
      self.f.__exit__(tp, value, traceback)


class GJob:
   """Class representing a gaussian job"""

   COMMANDPat = re.compile("\s*#[PTN]* ")
   SKIPGeomPat = re.compile("geom=\S*check", re.IGNORECASE)
   LINKPat  = re.compile("\s*--Link1--")
   OPTPat   = re.compile("\sopt([=\(]+([^\s\)]+)?\)?)", re.IGNORECASE)
   FREQPat  = re.compile("\sfreq[\s=\(]", re.IGNORECASE)
   MAXCylPat = re.compile("MaxCycles=\d+", re.IGNORECASE)
   CALCFCPat = re.compile("readFC|calcfc|calchffc|rcfc", re.IGNORECASE)
   GEOMPat   = re.compile("\s*geom(=\S+)", re.IGNORECASE)
   GUESSPat  = re.compile("\s*guess(=\S+)", re.IGNORECASE)

   def __init__(self, start, command, middle, coords, end):
      self.start = start       # % directives
      self.command = command   # gaussian command line
      self.middle = middle     # comment, charge and multiplicity
      self.coords = coords
      self.end = end           # anything after the coordinates

   def __str__(self):
      return ''.join(
         [self.start, self.command, self.middle, self.coords, self.end])

   def isOpt(self):
      return GJob.OPTPat.search(self.command)
   def isFreq(self):
      return GJob.FREQPat.search(self.command)

   def execute(self, outName):
      com = dedent("""
         date>>%s;gaussian.csh >>%s<<'gJOBComs'
         %s'gJOBComs'""") % (outName,outName, str(self))
      #warn(com)
      status = call(["/bin/csh", "-fc", com])
      if status > 0:
         raise IOError("Gaussian returned error code=%d" % status)
      p = Popen("tail -n 10 "+outName, shell=True, bufsize=2048, 
                stdin=PIPE, stdout=PIPE, close_fds=True)
      stdin,stdout= p.stdin, p.stdout
      #stdin,stdout = os.popen2("tail -n 10 "+outName)
      stdin.close()
      lines = stdout.read()
      stdout.close()

      return b" Normal termination of Gaussian" in lines

   def copy(self, chkGeom=False, optSteps='', optCalcFC=False, optReadFC=False):
      newCom = self.command
      newMiddle = self.middle
      newCoords = self.coords

      ma = GJob.OPTPat.search(newCom)
      if (optSteps or optCalcFC or optReadFC) and not ma:
         raise Exception("Not an optimization:" + str(self))
      elif optSteps or optCalcFC or optReadFC:
         optArgs= ma.group(2)
         if optSteps:
            optArgs= GJob.MAXCylPat.sub("",optArgs)
            if optArgs: optArgs += ","
            optArgs += "MaxCycles="+str(optSteps)
         if optCalcFC:
            optArgs = GJob.CALCFCPat.sub("",optArgs)
            if optArgs: optArgs += ","
            optArgs += "CalcFC"
         if optReadFC:
            optArgs = GJob.CALCFCPat.sub("",optArgs)
            if optArgs: optArgs += ","
            optArgs += "ReadFC"
         optArgs = optArgs.replace(",,",",")
         if optArgs.startswith(",") : optArgs = optArgs[1:]
         newCom = GJob.OPTPat.sub(" opt=(%s)"%optArgs,newCom)
      if chkGeom:
         newCom = GJob.GEOMPat.sub("",newCom)
         newCom = GJob.GUESSPat.sub("",newCom)
         newCom = newCom.rstrip() + " Geom=AllCheck Guess=TCheck\n"
         newMiddle = ""
         newCoords = ""
      return GJob(self.start, newCom, newMiddle, newCoords, self.end)

   @staticmethod
   def readNext(inFile):
      start = ""
      command = ""
      middle = ""
      coords = ""
      end = ""

      line = inFile.readline()
      if not line: return None
      while not GJob.COMMANDPat.match(line):
         start += line
         line = inFile.readRequiredLine()

      while line.strip():
         command += line
         line = inFile.readRequiredLine()

      if not GJob.SKIPGeomPat.search(command):
         middle = "\n"
         line = inFile.readRequiredLine()
         # read comment lines
         while line.strip():
            middle += line
            line = inFile.readRequiredLine()
         middle += line
         # read charge and multiplicity
         middle += inFile.readRequiredLine()

         line = inFile.readRequiredLine()
         while line.strip():
            coords += line
            line = inFile.readRequiredLine()

      while line and not GJob.LINKPat.match(line):
         end += line
         line = inFile.readline()

      return GJob(start, command, middle, coords, end)





desc = """Run guassian optimization run.
          Your gInFile may contain multiple jobs.
          Whenever an optimization job is found it will be executed in multiple
          subjobs with MaxCycle=optSteps. If the optimization does not
          complete a frequency calculation is done with the final geometry.
          If the n-1'ed step was a freq job it's parameters will be retained,
          if not then the "CalcFC" option will be added to the opt keyword.
Note that gOpt will modify your gaussian options somewhat.

Example: set n=myName.g ; set nCPU=4 ; mysub.py -q medium -jobName $n:r -nCPU $nCPU -totalMem 10 -- gOpt.py -in $n"""

parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
parser.add_argument('-in', dest='gInFileName', required=True,
       help='gaussian command file, out will be name.out')
parser.add_argument('-optSteps', dest='optSteps', required=False, default=8,
       help='Number of optimizaton steps to execute before recalculating freq (def=%d)'%8)
parser.add_argument('-restartJob', metavar="<n>", type=int, required=False, default=0,
       help='restart this computation with job number <n>. Only for opt jobs.')
args = parser.parse_args()

gInFileName = args.gInFileName
gOutFileName, dummy = os.path.splitext(gInFileName)
gOutFileName += ".out"
restartJob = args.restartJob
optSteps=args.__dict__.get('optSteps',8)

gJobs = []
with FileLineWrapper(open(gInFileName)) as gInFile:
   gJob = GJob.readNext(gInFile)
   while gJob:
      gJobs.append(gJob)
      gJob = GJob.readNext(gInFile)

lastGJob = None
for gJob in gJobs:

   restartJob -= 1
   if restartJob > 0:
       continue

   if gJob.isOpt(): # and lastGJob != None:
      newGJob = gJob.copy(optSteps=optSteps)
      success = newGJob.execute(gOutFileName)
      while not success:
         if lastGJob and lastGJob.isFreq():
            newGJob = lastGJob.copy(chkGeom=True)
            if not newGJob.execute(gOutFileName) :
               raise IOError("Freq calculation did not complete!")
            newGJob = gJob.copy(optSteps=optSteps,optReadFC=True)
            success = newGJob.execute(gOutFileName)
         else:
            newGJob = gJob.copy(chkGeom=True,optSteps=optSteps,optCalcFC=True)
            success = newGJob.execute(gOutFileName)
   else:
      gJob.execute(gOutFileName)
   lastGJob = gJob
