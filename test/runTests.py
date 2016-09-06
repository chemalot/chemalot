#!/usr/bin/env python
from __future__ import print_function, division
import xml.etree.ElementTree as ET
import glob
import sys
import re
import os
import os.path
import argparse
import subprocess
import shutil

desc= """run tests"""


def warn(*objs):
   print(*objs, file=sys.stderr)

dirRe = re.compile("\$dir")
outdirRe = re.compile("\$outDir")
tests = glob.glob("*/test.xml")

parser = argparse.ArgumentParser(description=desc)

parser.add_argument('-debug', dest='debug', action='store_true', default=False,
       help='Print additional information')

parser.add_argument(dest='tests', metavar='xmlPaths',
       default=tests, nargs='*',
       help='One or more full path to test.xml files to execute. '+
            'If none given all will be executed.')
args = parser.parse_args()
debug = args.debug

def replaceVars(dir, str):
   str = outdirRe.sub('$dir/out', str)
   str = dirRe.sub(dir, str)
   return str


def printResult(success, cpOpts, dir, wdir, out, err, ref):
   if not success:
      print("FAILED:\t%s %s" % (dir, out))
      if wdir :
         if not wdir.startswith("/"): 
            wdir = "%s/%s" % (dir,wdir)
         print("   wdir: %s" % wdir)
      else:
         wdir = dir
      if "-R" in cpOpts:
         if out == '.': out = ""  # out was just "." remove from path
         ref = ref + '/'
         out = out + '/*'  # copy all files from out to ref
      print("   stderr in: %s/%s" % (wdir,err))
      refPath = ref
      if not refPath.startswith("/") : refPath = wdir + "/" + ref
      print("   if correct issue: cp %s %s/%s %s" %
            (cpOpts, wdir, out, refPath))
   else:
      print("Success\t%s %s" % (dir, out))

startDir = os.getcwd()

for tXml in args.tests:
   if not tXml.endswith(".xml") : tXml = tXml + "/test.xml"
   dir=os.path.dirname(tXml)
   absDir=os.path.abspath(dir)
   shutil.rmtree(dir + '/out', True)
   os.mkdir(dir + '/out')
   try:
      root = ET.parse(tXml).getroot()

      # iterate over <test> elements
      for test in root.iter('test'):
         inp =replaceVars(absDir, test.attrib['in']  )
         out =replaceVars(absDir, test.attrib['out'] )
   
         # collect command from any innertext
         com = test.text.strip()
         for t in test.iter():
            com = com + " " + t.tail.strip()
         com =replaceVars(absDir, com.strip() )

         err=re.sub('\.[^.]+$','.err', out)

         if not inp or not out or not com:
            print("%s: inp, out or command missing! (%s, %s, %s)\n" 
                  % (test,inp,out,com))
            continue
      
         os.chdir(dir)

         wdir = ''
         if 'workDir' in test.attrib: 
            wdir = replaceVars( absDir, test.attrib['workDir'] )
            os.chdir(replaceVars(absDir, test.attrib['workDir']))

         errf = open( err, 'w+' )
   
         # runn all init scripts
         for init in test.iter('init'):
            icom = replaceVars(absDir, init.text.strip() )
            if debug: print("EXECUTING: %s" % icom)
            subprocess.call(["tcsh", "-fec", icom], stdout=errf, stderr=errf )

         # execute command
         # perl clears out 4th row in sdf file to avoid differences on timestamp
         com="cat %s \\\n|%s\\\n" % (inp, com)
         com="(%s)\\\n|perl -pe 's/^(  -OEChem-).+/$1/'" % (com)
         
         if debug: print("EXECUTING: %s" % com)
         sout = open( out, 'w+' )
         retCode = subprocess.call(["tcsh", "-fec", com],
                                   stdout=sout, stderr=errf)
         sout.close()
   
         if retCode != 0:
            print("FAILED with status %d:\t%s %s" % (retCode, dir, out))
            if wdir :
               if not wdir.startswith("/"): 
                  fwdir = "%s/%s" % (dir,wdir)
               print("   wdir: %s" % fwdir)
            else:
               fwdir = dir
            print("   stderr in: %s/%s" % (fwdir,err))
            os.chdir(startDir)
            continue
   
         # execute postprocess script
         for post in test.iter('postprocess'):
            pcom = replaceVars(absDir, post.text.strip() )
            if debug: print("EXECUTING: %s" % pcom)
            subprocess.call(["tcsh", "-fec", pcom], stdout=errf, stderr=errf )
   
         print("\n\n===========================================", file=errf)
         print("Diffs and checks", file=errf)
         print("===========================================", file=errf)
   
         # itereate over <diff> elements
         for diffs in test.iter('diff'):
            ref =replaceVars(absDir, diffs.attrib['ref'] )
            if debug: print("EXECUTING: diff %s %s" % (ref, out))
            ret = subprocess.call(["diff", ref, out], stdout=errf, stderr=errf )
            printResult( ret == 0, "", dir, wdir, out, err, ref)
   
         for diffs in test.iter('diffDir'):
            refDir = replaceVars(absDir, diffs.attrib['refDir'] )
            outDir = replaceVars(absDir, diffs.attrib['outDir'] )
            opts = diffs.attrib.get('opts', '') # opts is optional
            opts = replaceVars(absDir, opts )
            com = "diff %s '%s' '%s'" % (opts, refDir, outDir)
            if debug: print("EXECUTING: %s" % com)
            ret = subprocess.call(["tcsh", "-fec", com],
                                  stdout=errf, stderr=errf )
            printResult( ret == 0, "-R", dir, wdir, outDir, err, refDir);
   
         errf.close()
   
         os.chdir(startDir)
   except ET.ParseError as e:
      warn("\nError parsing: " + tXml)
      warn("   %s" % e)
