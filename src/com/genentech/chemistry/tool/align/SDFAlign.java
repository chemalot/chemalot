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

package com.genentech.chemistry.tool.align;

import java.util.ArrayList;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import openeye.oechem.OEAtomBase;
import openeye.oechem.OEAtomBaseIter;
import openeye.oechem.OEExprOpts;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEIFlavor;
import openeye.oechem.OEMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;
import openeye.oechem.oemolostream;

public class SDFAlign
{  public enum OUTType
   {  ALIGNED,ORIGINAL; }

   public static boolean is3DFormat(int fmt)
   {  if (fmt == OEFormat.SMI || fmt == OEFormat.ISM || fmt == OEFormat.CAN|| fmt == OEFormat.MF)
      {  return false;
      }
      return true;
   }

   public static void main(String...args)
   {  Options options = new Options();
      Option opt = new Option("in", true, "input sd file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out", true, "output file");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("method", true, "fss|sss|MCS|clique (default mcs).");
      options.addOption(opt);

      opt = new Option("ref", true, "reference molecule (if not given first molecule in '-in' is used). If multiple ref molecules are read the min RMSD is reported.");
      options.addOption(opt);

      opt = new Option("mirror", false, "If given and the molecule is not chiral, return best mirror image.");
      options.addOption(opt);

      opt = new Option("rmsdTag", true, "Tagname for output of rmsd, default: no output.");
      options.addOption(opt);
      
      opt = new Option("atomDeviationTag", true, "Tagname for output of csv with deviation of heavy atoms, default: no output.");
      options.addOption(opt);
      
      opt = new Option("atomMatch", true, "Sequence of none|default|hcount|noAromatic specifing how atoms are matched cf. oe document.\n"
                                         +"noAromatic can be used to make terminal atoms match aliphatic and aromatic atoms.\n"
                                         +"Queryfeatures are considered only if default is used.");
      options.addOption(opt);

      opt = new Option("bondMatch", true, "Sequence of none|default specifing how bonds are matched cf. oe document.");
      options.addOption(opt);

      opt = new Option("keepCoreHydrogens", false, "If not specified the hydrigen atoms are removed from the core.");
      options.addOption(opt);

      opt = new Option("outputMol", true, "aligned|original (def: aligned) use original to just compute rmsd.");
      options.addOption(opt);

      opt = new Option("doNotOptimize", false, "If specified the RMSD is computed without moving optimizing the overlay.");
      options.addOption(opt);

      opt = new Option("atomMapTags", true, "'|' separated list of field names in the ref file, "
                 + "containing ',' separated values per atom, "
                 + "After mapping the values are transfered to fields with the same name in the output file, "
                 + "reordered to be associated with the matched atoms. "
                 + "Note: only info on heavy atoms is transfered!");
      options.addOption(opt);

      opt = new Option("atomTypeTag", true, "tagName containing ',' separated list of atom types."
                 + " This is only used to verify that the atom types are in the same order as the atoms in the molecule."
                 + " Thus validating that the values in atomMapTags are ordered correctly." );
      options.addOption(opt);
      

      opt = new Option("quiet", false, "Reduced warining messages");
      options.addOption(opt);

      CommandLineParser parser = new BasicParser();
      CommandLine cmd = null;
      try
      {
         cmd = parser.parse(options, args);
      } catch (Exception e)
      {  exitWithHelp(e.getMessage(), options);
      }

      if (cmd.getArgs().length > 0)
         exitWithHelp("To many arguments", options);

      // do not check aromaticity on atoms so that a terminal atom matches aromatic and non aromatic atoms
      int atomExpr = OEExprOpts.DefaultAtoms;
      int bondExpr = OEExprOpts.DefaultBonds;

      String atomMatch = cmd.getOptionValue("atomMatch");
      if( atomMatch == null ) atomMatch = "";
      atomMatch = '|' + atomMatch.toLowerCase() + '|';

      String bondMatch = cmd.getOptionValue("bondMatch");
      if( bondMatch == null ) bondMatch = "";
      bondMatch = '|' + bondMatch.toLowerCase() + '|';

      String inFile = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String refFile = cmd.getOptionValue("ref");
      String method = cmd.getOptionValue("method");
      String atomDevTag = cmd.getOptionValue("atomDeviationTag");
      String rmsdTag = cmd.getOptionValue("rmsdTag");
      String oMol = cmd.getOptionValue("outputMol");
      boolean doMirror = cmd.hasOption("mirror");
      boolean doOptimize = !cmd.hasOption("doNotOptimize");
      boolean quiet = cmd.hasOption("quiet");

      OUTType outputMol = oMol == null ? OUTType.ALIGNED : OUTType.valueOf(oMol.toUpperCase());
      
      if( atomDevTag != null && method != null && method.equals("fss") )
            throw new Error("-atomDeviationTag is not supported for method = fss");
      
      String[] atomInfoTags = new String[0];
      if( cmd.hasOption("atomMapTags"))
         atomInfoTags = cmd.getOptionValue("atomMapTags").split("\\|");
      if( atomInfoTags.length > 0 && method != null && method.equals("fss") )
         throw new Error("-atomMapTags is not supported for method = fss");
      if( refFile == null )
         throw new Error("-atomMapTags is only works with -ref");
      String atomTypeTag = cmd.getOptionValue("atomTypeTag"); 
   
      if( atomMatch.startsWith("|none") )      atomExpr = 0;
      if( atomMatch.contains("|hcount|") )     atomExpr |= OEExprOpts.HCount;
      if( atomMatch.contains("|noAromatic|") ) atomExpr &= (~ OEExprOpts.Aromaticity);

      if( bondMatch.startsWith("|none") ) bondExpr = 0;

      ArrayList<OEMol> refmols = getReferences(cmd.hasOption("keepCoreHydrogens"), refFile, atomInfoTags, atomTypeTag);
            
      oemolistream fitfs = new oemolistream(inFile);
      if (!is3DFormat(fitfs.GetFormat()))
         throw new Error("Invalid input format: need 3D coordinates");

      oemolostream ofs = new oemolostream(outFile);
      if (!is3DFormat(ofs.GetFormat()))
         throw new Error("Invalid output format: need 3D coordinates");

      AlignInterface aligner = null;
      OEGraphMol fitmol = new OEGraphMol();

      if( oechem.OEReadMolecule(fitfs, fitmol) )
      {  if( refmols.size() == 0)
         {  OEMol rmol = new OEMol(fitmol);
            if( ! cmd.hasOption("keepCoreHydrogens"))
               oechem.OESuppressHydrogens(rmol);

            refmols.add(rmol);
         }

         if( "sss".equalsIgnoreCase(method))
         {  aligner = new SSSAlign(refmols, outputMol, rmsdTag, atomDevTag, atomInfoTags, doOptimize, doMirror, atomExpr, bondExpr, quiet);

         } else if( "clique".equalsIgnoreCase(method))
         {  aligner = new CliqueAlign(refmols, outputMol, rmsdTag, atomDevTag, atomInfoTags, doOptimize, doMirror, atomExpr, bondExpr, quiet);

         }else if( "fss".equalsIgnoreCase(method))
         {  if( cmd.hasOption("atomMatch") || cmd.hasOption("bondMatch") )
               exitWithHelp("method fss does not support '-atomMatch' or '-bondMatch'", options);
            aligner = new FSSAlign(refmols, outputMol, rmsdTag, doOptimize, doMirror);

         } else
         {  aligner = new McsAlign(refmols, outputMol, rmsdTag, atomDevTag, atomInfoTags, doOptimize, doMirror, atomExpr, bondExpr, quiet);
         }

         do
         {  aligner.align(fitmol);
            oechem.OEWriteMolecule(ofs, fitmol);
         }while (oechem.OEReadMolecule(fitfs, fitmol));

      }

      fitmol.delete();
      if( aligner != null ) aligner.close();
      for( OEMolBase mol : refmols) mol.delete();
      fitfs.close();
      ofs.close();
   }

   protected static ArrayList<OEMol> getReferences(boolean keepCoreHydrogens, String refFile, 
                                                   String[] atomInfoTags, String atomTypeTag) throws Error
   {  ArrayList<OEMol> refmols = new ArrayList<OEMol>();
    
      if( refFile != null )
      {  oemolistream reffs = new oemolistream(refFile);
         if (!is3DFormat(reffs.GetFormat()))
            throw new Error("Invalid input format: need 3D coordinates");
         reffs.SetFormat(OEFormat.MDL);

         int aromodel = OEIFlavor.Generic.OEAroModelOpenEye;
         int qflavor  = reffs.GetFlavor( reffs.GetFormat());
         reffs.SetFlavor(reffs.GetFormat(),(qflavor|aromodel));

         OEMol rmol = new OEMol();
         while(oechem.OEReadMDLQueryFile(reffs,rmol))
         {  attachAtomData(rmol, atomInfoTags, atomTypeTag);
         
            if( ! keepCoreHydrogens )
               oechem.OESuppressHydrogens(rmol);
            refmols.add(rmol);
            rmol = new OEMol();
         }
         rmol.delete();

         if( refmols.size() == 0)
            throw new Error("reference file had no entries");

         reffs.close();
      }
      return refmols;
   }

   /**
    * the tags in atomInfoTags contain "," separated values in atom order. Add them 
    * as StringData to the atoms so they can be transfered to the matched atoms after
    * alignment.
    * 
    * @param rmol refernce molecules
    * @param atomInfoTags the tags point to "," separated values in atom order
    * @param atomTypeTag field name of atom label in atom order. 
    *                    Used to validate that atom and data order was not mixed up.
    */
   private static void attachAtomData(OEMol rmol, String[] atomInfoTags, String atomTypeTag)
   {  checkAtomTypeTagSequence(rmol, atomTypeTag);
   
      for( String tag : atomInfoTags)
      {  String atomTagVal = oechem.OEGetSDData(rmol, tag);
         if( atomTagVal != null && atomTagVal.length() > 0)
         {  String[] atomTagVals = atomTagVal.split(" *, *",-1);
            if( atomTagVals.length != rmol.NumAtoms()) 
               System.err.print(
                  String.format("mol '%s' tag=%s number of values does not match atom count!\n",
                     rmol.GetTitle(), tag ));
            
            OEAtomBaseIter atit = rmol.GetAtoms();
            while( atit.hasNext())
            {  OEAtomBase at = atit.next();
               int idx = at.GetIdx();
               if( idx < atomTagVals.length)
                  at.SetStringData(tag, atomTagVals[idx]);
            }
            atit.delete();
         }
      }
   }

   /**
    * Ensure that atom types in comma separated field atomTypeTag are in same order as 
    * the atom types of this molecule.
    * 
    * This is a check to alert us if atoms have been reordered as compared to the
    * values stored in the fields from  atomInfoTags
    */
   private static void checkAtomTypeTagSequence(OEMol rmol, String atomTypeTag) throws Error
   {
      if( atomTypeTag != null && atomTypeTag.length() > 0 )
      {  String atomTypesCSV = oechem.OEGetSDData(rmol, atomTypeTag);
         if( atomTypesCSV == null || atomTypesCSV.length() == 0)
            return;
         
         String[] atomTypes = atomTypesCSV.toLowerCase().split(" *, *",-1);
      
         if( rmol.NumAtoms() < atomTypes.length )
            throw new Error(String.format("molecule (%s) has less atoms (%d) than listed in atomTypes (%s)",
                     rmol.GetTitle(), rmol.NumAtoms(), atomTypesCSV));
         
         if( rmol.NumAtoms() > atomTypes.length )
            System.err.printf("molecule (%s) has more atoms (%d) than listed in atomTypes (%s)\n",
                     rmol.GetTitle(), rmol.NumAtoms(), atomTypesCSV);
               
         OEAtomBaseIter atIt = rmol.GetAtoms();
         while( atIt.hasNext() )
         {  OEAtomBase at = atIt.next();
            if(at.GetIdx() >= atomTypes.length ) continue;
            
            String sym = oechem.OEGetAtomicSymbol(at.GetAtomicNum()).toLowerCase();
            if(atomTypes[at.GetIdx()].length() > 0 && ! sym.equals( atomTypes[at.GetIdx()]))
            {  System.err.printf("molecule (%s) has more mismatch in atoms listed in %s, atom (%d) '%s'!='%s'\n",
                  rmol.GetTitle(), atomTypeTag, at.GetIdx(), sym, atomTypes[at.GetIdx()]);
               atIt.delete();
               return;
            }
         }
         atIt.delete();
      }
   }

   private static void exitWithHelp(String msg, Options options)
   {  System.err.println(msg);

      HelpFormatter formatter = new HelpFormatter();
      formatter.setWidth(120);
      formatter.printHelp( "sdfAlign -in .sdf -out .sdf [-ref .sdf]", options );
      System.exit(1);
   }

}
