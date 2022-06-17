/*
   Copyright 2008-2015 Genentech Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/
package com.genentech.chemistry.openEye.apps;

import java.io.*;
import java.util.*;

import openeye.oechem.*;

import org.apache.commons.cli.*;

import com.aestel.io.IOUtil;
import com.genentech.oechem.tools.OETools;

/**
 * Create gaussian (in the future other) input for QM torsion scan.
 * Define torsion by reading file with 4 atoms. The two central atoms define the
 * bond around which to scan the torsion.
 *
 * Then read sdf file and find the 4 atoms closest to torsion atoms. Create n
 * conformations which differ in torsion angle and write one gaussian input file for each.
 *
 * @author Man-Ling Lee / August 04, 2012 Copyright 2012 Genentech
 */
public class QTorsionProfileGenerator
{  public static final int FRAGNumTag = oechem.OEGetTag("fragNum");


   private static final String MY_NAME = QTorsionProfileGenerator.class.getSimpleName();
   private static final String OPT_INFILE = "in";
   private static final String OPT_SDFFILE = "sdf";
   private static final String OPT_OUTPREFIX = "outPrefix";
   private static final String OPT_TEMPLATE = "template";
   private static final String OPT_BONDFILE = "bondFile";
   private static final String OPT_STARTTorsion = "startTorsion";
   private static final String OPT_TORSIONIncrement = "torsionIncrement";
   private static final String OPT_NSTEPS = "nSteps";
   private static final String OPT_COREFILE = "core";
   private static final String OPT_OUTNAMETAG = "outNameTag";
   private static final String OPT_WORKDIR = "workDir";
   private static final String OPT_NCPU = "nCPU";
   private static final String OPT_MINIMIZE = "minimize";
   private static final String OPT_CONSTRIANT = "constraintStrength";
   private static final String OPT_MAXCONFS_PER_STEP = "maxConfsPerStep";
   private static final String OPT_MEM = "mem";
   private static final String OPT_DEBUG = "debug";
   private static final String OPT_SPIN_MULTIPLICITY = "spinMultiplicity";



   private static final String FROZEN_ATOM_TAG = "qTorFrozen";


   private final String outPrefix;
   private final double startTorsion;
   private final double torsionIncrement;
   private final int nSteps;
   private final boolean doMinimize;
   private final String qTemplate;
   private final String outNameTag;
   private final String workDir;
   private final String bondFile;
   private final int  sampleNConf;
   private final String constraintStrength;
   private final boolean debug;
   private final String spinMult;


   /**
    *
    * @param templateName template file name for QM input, contains
    *                      #XYZ#, #FNAME# and #NCPU# place holders
    * @param outPrefix prefix for QM output files must be null if outNameTag not null
    * @param outNameTag tag name for sdf tag containing prefix for QM input file
    * @param bondFile file containing the 4 atoms defining the torsion.
    * @param startTorsion start angle for torsion, if NaN use torsion as in input
    * @param torsionIncrements increment in deg for each step
    * @param nSteps number of steps = number of QM command files per sdf record
    * @param doMinimize
    * @param constraintStrength one of strong (90),medium (45), weak(20), none or a floating point number specifying the strength of tethered constrains for -doMinimize
    */
   private QTorsionProfileGenerator(String templateName,
            String workDir, String outPrefix, String outNameTag,
            String bondFile, String spinMult, double startTorsion,
            double torsionIncrements, int nSteps, int sampleNConf,
            boolean doMinimize, String constraintStrength, boolean debug)
   {  this.workDir = workDir;
      this.outPrefix = outPrefix;
      this.outNameTag = outNameTag;
      this.startTorsion = startTorsion;
      this.torsionIncrement = torsionIncrements;
      this.nSteps = nSteps;
      this.doMinimize = doMinimize;
      this.bondFile = bondFile;
      this.spinMult = spinMult;
      this.sampleNConf = sampleNConf;
      this.constraintStrength = constraintStrength;
      this.debug = debug;

      String dummy;
      try
      {  dummy = IOUtil.fileToString(templateName);
      } catch (IOException e)
      {  throw new Error("Template file could not be read");
      }

      // add "modredundant" option to "opt" keyword in gaussian command
      dummy = dummy.replaceAll("(#.*\\bopt\\b=?\\(?)((?!\\S*ModRedundant)\\S+)", "$1ModRedundant,$2");
      this.qTemplate = dummy;
      if(! qTemplate.contains("#XYZ#") )
         throw new Error("Template file does not contain placeholder #XYZ#");
   }



   /**
    * Create gaussian files and core file for one sdf file with molecules.
    *
    * @param inFile name of input file with molecules
    * @param coreFile name of file into which to store the core, may be null
    * @param nCPU overwrite gaussian %nprocshared parameter, if null no overwriting takes place
    * @param sdfFile if null no sdf file is written
    * @param memStr Memory for Gaussian e.g. "10GB"
    */
   private void run(String inFile, String sdfFile, String coreFile, String nCPU, String memStr) throws IOException
   {  TorsionScanner torScanner = new TorsionScanner(bondFile, coreFile, FROZEN_ATOM_TAG,
                                      nSteps, startTorsion, torsionIncrement,
                                      sampleNConf, doMinimize, constraintStrength);


      oemolithread ifs = new oemolithread(inFile);
      OEGraphMol mol = new OEGraphMol();
      List<OEGraphMol> inMols =  new ArrayList<OEGraphMol>();

      oemolothread ofs = null;
      if( sdfFile != null ) ofs = new oemolothread(sdfFile);

      int molNum = 0;
      while (oechem.OEReadMolecule(ifs, mol))
      {  try
         {  Map<String,Integer> angleCountMap = new HashMap<String,Integer>();

            if(coreFile != null) inMols.add(new OEGraphMol(mol));
            OEMCMolBase mcMol = torScanner.run(mol);
            OEConfBaseIter cIt = mcMol.GetConfs();
            while( cIt.hasNext() )
            {  OEConfBase conf = cIt.next();
               if( ofs != null ) oechem.OEWriteMolecule(ofs, conf);
               String angle = oechem.OEGetSDData(conf, TorsionScanner.ANGLE_TAG+"_1");
               String angle2 = oechem.OEGetSDData(conf, TorsionScanner.ANGLE_TAG + "_2");
               
               String sameAngleStr = angle;
               if( angle2 != null && angle2.length() > 0)
               {  sameAngleStr += '_' + angle;
                  angle = String.format("%04.0f_%04.0f", Double.parseDouble(angle), Double.parseDouble(angle2));
               }else
               {  angle = String.format("%04.0f", Double.parseDouble(angle));
               }
               
               int countSameAngle = 0;
               if( angleCountMap.containsKey(sameAngleStr) )
                  countSameAngle = angleCountMap.get(sameAngleStr);
               angleCountMap.put(sameAngleStr,  ++countSameAngle);

               molNum++;

               String fName;
               if( outPrefix != null )
                  fName = String.format("%s_%03d.%s_%03d", outPrefix, molNum, angle, countSameAngle);
               else if( "TITLE".equals(outNameTag))
                  fName = String.format("%s.%s_%03d", mol.GetTitle(), angle, countSameAngle);
               else
                  fName = String.format("%s.%s_%03d", oechem.OEGetSDData(mol, outNameTag), angle, countSameAngle);

               String xMat = getXMat(conf);

               String frozenStmt = getFrozenCoordinate(conf);
               String gCom = qTemplate.replace("#FName#", fName);

               if( qTemplate.contains("#FIX#") )
                  gCom = gCom.replace("#FIX#", frozenStmt);
               else
                  // add statement to freeze torsion
                  xMat += "\n" + getFrozenCoordinate(conf);

               gCom = gCom.replace("#XYZ#", xMat);
               gCom = gCom.replace("#mem#", memStr);
               if( nCPU != null )
                  gCom = gCom.replaceAll("%nprocshared=.*", "%nprocshared="+nCPU);

               IOUtil.stringToFile(workDir + File.separatorChar + fName+".g", gCom);
            }
            cIt.delete();
            mcMol.delete();

         } catch(IllegalArgumentException e)
         {  System.err.println(e.getMessage() + " compound ignored!\n"+ e.getMessage());
            if( debug ) e.printStackTrace(System.err);
         }
      }
      mol.delete();

      if( ofs != null )
      {  ofs.close();
         ofs.delete();
      }

      ifs.close();
      ifs.delete();

      if(coreFile != null && molNum > 0)
         torScanner.computeCore(coreFile);

      torScanner.close();
   }


   private static String getFrozenCoordinate(OEConfBase conf)
   {  String res= getFrozenCoordinate(conf, FROZEN_ATOM_TAG);
      
      if( oechem.OEHasSDData(conf, FROZEN_ATOM_TAG + "_2"))
         res = res + getFrozenCoordinate(conf, FROZEN_ATOM_TAG + "_2");

      return res;
   }

   private static String getFrozenCoordinate(OEConfBase conf, String frozenAtomTag)
   {  String[] fAtomStr = oechem.OEGetSDData(conf, frozenAtomTag).split(" ");

      // increment atoms by one because guassian uses 1 based indexes
      StringBuilder sb = new StringBuilder(fAtomStr.length+4);
      for(String s: fAtomStr)
         sb.append(Integer.parseInt(s)+1).append(" ");

      return String.format("D %s F\n", sb.toString());
   }


   private String getXMat(OEMolBase mol)
   {
      StringBuilder sb = new StringBuilder();
      sb.append(OETools.getCharge(mol)).append(' ').append(spinMult).append("\n");

      float coords[] = new float[mol.GetMaxAtomIdx() * 3];
      mol.GetCoords(coords);

      OEAtomBaseIter atIt = mol.GetAtoms();
      while (atIt.hasNext())
      {
         OEAtomBase at = atIt.next();
         int idx = at.GetIdx();
         sb.append(String.format(" %s     %.4f %.4f %.4f\n",
                  oechem.OEGetAtomicSymbol(at.GetAtomicNum()), coords[idx * 3],
                  coords[idx * 3 + 1], coords[idx * 3 + 2]));
      }
      atIt.delete();
      return sb.toString();
   }




   public static void main(String... args) throws IOException
   { // create command line Options object
      Options options = new Options();
      Option opt = new Option(OPT_INFILE, true,
               "input file oe-supported Use .sdf|.smi to specify the file type.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_SDFFILE, true,
               "file to write onformers for debugging (oe-supported Use .sdf|.smi to specify the file type).");
      options.addOption(opt);

      opt = new Option(OPT_WORKDIR, true, "Write files into this directory!");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_OUTPREFIX, true,
               "Prefix for output file.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_OUTNAMETAG, true,
               "TagName of field containing outFilePrefix. (Use TITLE for mol title).");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_TEMPLATE, true,
               "Template file for quantum program containing #XYZ# place holder line. A #FName# placeholder can be used fro chk files and the like.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_SPIN_MULTIPLICITY, true, "Spin Multiplicity of the input molecules (default 1)");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_NCPU, true, "Overwrite nprocshared parameter in guassian input file if given");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_MEM, true, "Memory for gaussian default='10GB' replaces #mem# in tempalte");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_MINIMIZE, false, "minimize conformer at each step using MMFFs");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_CONSTRIANT, true, "one of strong (90),medium (45), weak(20), none or a floating point number"
                                            +" specifying the strength of tethered constrains for -doMinimize (def=strong)");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(
               OPT_BONDFILE,
               true,
               "Structure file containing 4-8 atoms defining 1 oer 2 torsions. The first bond is a2-a3 the second is a-2-a-3."
              +"In each input molecule the atoms colses these atoms are used to define the torsion.");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_STARTTorsion, true,
               "The torsion in your inMol will be rotated by this value for the first job");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_TORSIONIncrement, true,
               "Incremnt each subsequent conformation by this step size");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_NSTEPS, true, "Number of conformations to create");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option(OPT_MAXCONFS_PER_STEP, true, "While holding the torsion fixed, maximum number of conformations of free atoms to generate.  default=1");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_COREFILE, true, "Outputfile to store guessed core.");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option(OPT_DEBUG, false, "Produce more debug output.");
      opt.setRequired(false);
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse(options, args);
      } catch (Exception e)
      {
         System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      args = cmd.getArgs();
      if (args.length != 0)
      {  System.err.println("Unknown arguments" + args);
         exitWithHelp(options);
      }

      if (cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      String inFile    = cmd.getOptionValue(OPT_INFILE);
      String sdfFile   = cmd.getOptionValue(OPT_SDFFILE);
      String outPrefix = cmd.getOptionValue(OPT_OUTPREFIX);
      String outNameTag= cmd.getOptionValue(OPT_OUTNAMETAG);
      String workDir   = cmd.getOptionValue(OPT_WORKDIR);
      String coreFile  = cmd.getOptionValue(OPT_COREFILE);
      String tempalte  = cmd.getOptionValue(OPT_TEMPLATE);
      String bondFile  = cmd.getOptionValue(OPT_BONDFILE);
      String nCPU      = cmd.getOptionValue(OPT_NCPU);
      String memStr    = cmd.getOptionValue(OPT_MEM);
      String maxConfsStr=cmd.getOptionValue(OPT_MAXCONFS_PER_STEP);
      String spinMult   =cmd.getOptionValue(OPT_SPIN_MULTIPLICITY);


      boolean doMinimize=cmd.hasOption(OPT_MINIMIZE);
      String constraintStrength = cmd.getOptionValue(OPT_CONSTRIANT);
      int nStep = Integer.parseInt(cmd.getOptionValue(OPT_NSTEPS));

      if(memStr == null || memStr.length() == 0 ) memStr = "10GB";

      if( spinMult == null || spinMult.length() == 0 ) spinMult = "1";

      if((outPrefix == null && outNameTag == null) || (outPrefix != null && outNameTag != null))
      {  System.err.println("Exactly one of -outPrefix or outNameTag must be given!");
         exitWithHelp(options);
      }

      if( workDir == null || workDir.trim().length() == 0 ) workDir = ".";

      double startTorsion = Double.NaN;
      if( cmd.hasOption(OPT_STARTTorsion))
         startTorsion = Double.parseDouble(cmd.getOptionValue(OPT_STARTTorsion));

      double torInc = Double.parseDouble(cmd.getOptionValue(OPT_TORSIONIncrement));
      int maxStepConfs = maxConfsStr == null ? 1 : Integer.parseInt(maxConfsStr);


      QTorsionProfileGenerator calculator = new QTorsionProfileGenerator(
               tempalte, workDir, outPrefix, outNameTag, bondFile, spinMult,
               startTorsion, torInc, nStep, maxStepConfs, doMinimize, constraintStrength,
               cmd.hasOption(OPT_DEBUG));

      calculator.run(inFile, sdfFile, coreFile, nCPU, memStr);
   }

   private static void exitWithHelp(Options options)
   {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp(MY_NAME, options);
      System.exit(1);
   }

}
