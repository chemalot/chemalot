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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import openeye.oechem.*;

import org.apache.commons.cli.*;

/**
 * SphereExclusion method implemented using internal fp toolkit and OEchem sdf reader.
 *
 * @author albertgo
 *
 */
public class SDFSubRMSD
{
   private SDFSubRMSD()
   {
   }


   public static void main(String...args) throws IOException
   {  // create command line Options object
      Options options = new Options();
      Option opt = new Option("in",true, "input file oe-supported");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("fragFile",true, "file with single 3d substructure query");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("isMDL",false, "if given the fragFile is suposed to be an mdl query file, query features are supported.");
      opt.setRequired(false);
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp(options);
      }
      args = cmd.getArgs();

      if(cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }

      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      String fragFile = cmd.getOptionValue("fragFile");

      // read fragment
      OESubSearch ss;
      oemolistream ifs = new oemolistream(fragFile);
      OEMolBase mol;
      if( ! cmd.hasOption("isMDL"))
      {  mol = new OEGraphMol();
         oechem.OEReadMolecule(ifs, mol);
         ss = new OESubSearch(mol, OEExprOpts.AtomicNumber, OEExprOpts.BondOrder);
      }
      else
      {  int aromodel = OEIFlavor.Generic.OEAroModelOpenEye;
         int qflavor  = ifs.GetFlavor(ifs.GetFormat());
         ifs.SetFlavor(ifs.GetFormat(),(qflavor|aromodel));
         int opts = OEMDLQueryOpts.Default|OEMDLQueryOpts.SuppressExplicitH;
         OEQMol qmol = new OEQMol();
         oechem.OEReadMDLQueryFile(ifs,qmol,opts);
         ss = new OESubSearch(qmol);
         mol = qmol;
      }

      double nSSatoms = mol.NumAtoms();
      double sssCoords[] = new double[mol.GetMaxAtomIdx() * 3];
      mol.GetCoords(sssCoords);
      mol.Clear();
      ifs.close();

      if(! ss.IsValid())
         throw new Error("Invalid query " + args[0]);


      ifs = new oemolistream(inFile);
      oemolostream ofs = new oemolostream(outFile);
      int count = 0;

      while (oechem.OEReadMolecule(ifs, mol))
      {  count++;
         double rmsd = Double.MAX_VALUE;
         double molCoords[] = new double[mol.GetMaxAtomIdx() * 3];
         mol.GetCoords(molCoords);

         for (OEMatchBase mb : ss.Match(mol, false))
         {  double r = 0;
            for (OEMatchPairAtom mp : mb.GetAtoms())
            {  OEAtomBase asss = mp.getPattern();
               double sx = sssCoords[asss.GetIdx()*3];
               double sy = sssCoords[asss.GetIdx()*3];
               double sz = sssCoords[asss.GetIdx()*3];

               OEAtomBase amol = mp.getTarget();
               double mx = molCoords[amol.GetIdx()*3];
               double my = molCoords[amol.GetIdx()*3];
               double mz = molCoords[amol.GetIdx()*3];

               r += Math.sqrt( (sx-mx)*(sx-mx) + (sy-my)*(sy-my) + (sz-mz)*(sz-mz) );
            }
            r /= nSSatoms;
            rmsd = Math.min(rmsd, r);
         }

         if( rmsd != Double.MAX_VALUE )
            oechem.OESetSDData(mol, "SSSrmsd", String.format("%.3f", rmsd));

         oechem.OEWriteMolecule(ofs, mol);
         mol.Clear();
      }

      ifs.close();
      ofs.close();

      mol.delete();
      ss.delete();
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( "SDFSubRMSD adds sub rmsd of substructure match to input structure."
                          +" Inputs must have been aligned to sss", options );
      System.exit(1);
   }
}

