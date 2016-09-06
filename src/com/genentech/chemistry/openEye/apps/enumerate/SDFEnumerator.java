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
package com.genentech.chemistry.openEye.apps.enumerate;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import openeye.oechem.*;
import openeye.oedepict.oedepict;

import org.apache.commons.cli.*;

import com.aestel.io.IOUtil;
import com.genentech.oechem.tools.OETools;

/**
 *
   -out c:\tmp\enum.smi
   ([Ag][*:2].[Be][*:3].[Ce][*:5]).([Ag][*:1]).([Be][*:4]).([Ce][*:6])>>[*:1]-[*:2].[*:3]-[*:4].[*:5]-[*:6]
   code\java\src\com\genentech\chemistry\openEye\apps\enumerate\core.smi
   code\java\src\com\genentech\chemistry\openEye\apps\enumerate\R1.smi
   code\java\src\com\genentech\chemistry\openEye\apps\enumerate\R2.smi
   code\java\src\com\genentech\chemistry\openEye\apps\enumerate\R3.smi
 */
public class SDFEnumerator
{  private final int nPositions;
   private final OELibraryGen libgen;
   private final LinkedHashMap<String, OEMolBase>[] reactants;
   private int counter=1;
   private final String uniqueTag = "sdfEnumerator_counter";
   private final String enumerated = "Enumerated";

   @SuppressWarnings("unchecked")
   SDFEnumerator(OELibraryGen libgen, boolean reactAllSites, String ... smiOrFiles)
   {  this.libgen = libgen;
      this.nPositions = smiOrFiles.length;
      this.reactants = new LinkedHashMap[nPositions];
      for( int i=0; i< nPositions; i++)
      {  readSmiOrFile( smiOrFiles[i], i, reactAllSites );
      }

   }


   private void readSmiOrFile(String smiOrFile, int pos, boolean reactAllSites)
   {  reactants[pos] = new LinkedHashMap<String, OEMolBase>();

      OEMolBase mol = new OEGraphMol();

      // check to see if smiOrFile is a valid smiles
      if( ! smiOrFile.startsWith(".") && ! new File(smiOrFile).exists() )
      {  // check if we got a smiles instead of a filename
         OETools.smiToMol(mol, smiOrFile);
         if( mol.IsValid() && mol.NumAtoms() > 0)
         {  oechem.OEAddSDData(mol, uniqueTag+pos, Integer.toString(counter));
            reactants[pos].put(Integer.toString(counter), new OEGraphMol(mol)); // create a copy to be outputted later
            counter++;
            libgen.AddStartingMaterial(mol, pos, ! reactAllSites );
            mol.delete();

            return;
         }
         mol.Clear();
      }

      // smiOrFile must be a filename
      oemolithread ifs = new oemolithread(smiOrFile);
      while(oechem.OEReadMolecule(ifs, mol))
      {  oechem.OEAddSDData(mol, uniqueTag+pos, Integer.toString(counter));
         reactants[pos].put(Integer.toString(counter), new OEGraphMol(mol)); // create a copy to be outputted later
         counter++;
         libgen.AddStartingMaterial(mol, pos, ! reactAllSites );
      }
      mol.delete();
      ifs.close();
      ifs.delete();
   }

   public void generateLibrary(final String outName,
                                 final int maxAtoms, final double randomFract, final boolean regenerate2D, final String unreactedFile)
   {  oemolothread ofs = new oemolothread(outName);
      for (OEMolBase mol : libgen.GetProducts())
      {  if( maxAtoms > 0 ) oechem.OESuppressHydrogens(mol);

         if (regenerate2D == true) oedepict.OEPrepareDepiction(mol, true, true);

         if( (randomFract > 1D || randomFract >= Math.random())
          && (maxAtoms == 0 || mol.NumAtoms() <= maxAtoms) )
            oechem.OEAddSDData(mol, enumerated, "Yes");
            oechem.OEWriteMolecule(ofs, mol);

            if (unreactedFile != null)
            {  // get uniqueTag of mol, remove that entry from reactants hashmap

               for(int i=0; i<nPositions; i++)
               {  String mol_key = oechem.OEGetSDData(mol, uniqueTag+i);
                  OEMolBase r = reactants[i].get(mol_key);
                  if (r != null)
                  {  reactants[i].remove(mol_key);
                     r.delete();
                  }
               }
            }
         mol.delete();
      }

      //print unreacted reagents
      if (unreactedFile != null)
      {  oemolothread ofsUnreacted = new oemolothread(unreactedFile);
         for(int i=0; i<nPositions; i++)
         {  for (OEMolBase mol : reactants[i].values())
            {  oechem.OEAddSDData(mol, enumerated, "No");
               oechem.OEWriteMolecule(ofsUnreacted, mol);
               mol.delete();
            }
         }
         ofsUnreacted.close();
         ofsUnreacted.delete();
      }

      ofs.close();
      ofs.delete();
   }

   public static void main(String ...args) throws IOException
   {  Options options = new Options();
      Option opt = new Option("out",true, "output file oe-supported");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("hydrogenExplicit",false, "Use explicit hydrogens");
      options.addOption(opt);

      opt = new Option("correctValences",false, "Correct valences after the enumeration");
      options.addOption(opt);

      opt = new Option("regenerate2D",false, "Regenerate 2D coordinates for the products");
      options.addOption(opt);

      opt = new Option("reactAllSites",false, "Generate a product for each match in a reagent.");
      options.addOption(opt);

      opt = new Option("randomFraction",true, "Only output a fraction of the products.");
      options.addOption(opt);

      opt = new Option("maxAtoms",true, "Only output products with <= maxAtoms.");
      options.addOption(opt);

      opt = new Option("notReacted", true, "Output file for reagents that didn't produce at leaste one output molecule, useful for debugging SMIRKS.");
      options.addOption(opt);

      CommandLineParser parser = new PosixParser();
      CommandLine cmd = null;
      try
      {  cmd = parser.parse( options, args);
      } catch(Exception e)
      {  System.err.println(e.getMessage());
         exitWithHelp("", options);
      }

      args = cmd.getArgs();
      if( args.length < 2 )
      {  exitWithHelp("Transformation and/or reagentFiles missing", options);
      }
      String smirks = args[0];
      if( new File(smirks).canRead() )
         smirks = IOUtil.fileToString(smirks).trim();
      if( ! smirks.contains(">>") )
         smirks = scaffoldToSmirks(smirks);
      String[] reagentSmiOrFiles = Arrays.copyOfRange(args, 1, args.length);

      if(cmd.hasOption("d"))
      {  System.err.println("Start debugger and press return:");
         new BufferedReader(new InputStreamReader(System.in)).readLine();
      }


      String outFile = cmd.getOptionValue("out");
      OELibraryGen lg = new OELibraryGen();
      lg.Init(smirks);
      if( ! lg.IsValid() )
         exitWithHelp("Invalid Transform: " + smirks, options);

      lg.SetExplicitHydrogens(cmd.hasOption("hydrogenExplicit"));
      lg.SetValenceCorrection(cmd.hasOption("correctValences"));
      lg.SetRemoveUnmappedFragments(true);

      boolean regenerate2D = cmd.hasOption("regenerate2D");
	   boolean reactAllSites = cmd.hasOption("reactAllSites");
      String unreactedFile=null;
	   if (cmd.hasOption("notReacted")) {
	      unreactedFile = cmd.getOptionValue("notReacted");
       }

      double randomFract = 2;
      if( cmd.hasOption("randomFraction" ) )
         randomFract = Double.parseDouble(cmd.getOptionValue("randomFraction"));

      int maxAtoms = 0;
      if( cmd.hasOption( "maxAtoms" ))
         maxAtoms = Integer.parseInt(cmd.getOptionValue("maxAtoms"));

      SDFEnumerator en = new SDFEnumerator(lg, reactAllSites, reagentSmiOrFiles);
      en.generateLibrary( outFile, maxAtoms, randomFract, regenerate2D, unreactedFile);
      en.delete();
   }

   private static final int SMIFlags = OESMILESFlag.AtomStereo|OESMILESFlag.BondStereo
            |OESMILESFlag.AtomMaps  |OESMILESFlag.Isotopes;


   private static String scaffoldToSmirks(String scaffoldSmi)
   {  // Goal convert [U+1]c1nc([U+2])ncc1 to
      // [U+101][c:1]1[n:2][c:3]([U+102])[n:4][c:5][c:6]1.[U+][*:7].[U+2][*:8];
      // >> [*:7][c:1]1[n:2][c:3]([*:8])[n:4][c:5][c:6]1
      OEMolBase mol = new OEGraphMol();
      OETools.smiToMol(mol, scaffoldSmi);
      if( ! mol.IsValid() || mol.NumAtoms() == 0 )
         System.err.println("Invalid smiles: " + scaffoldSmi);

      // make sure rgroups are ordered by fgroup number == charge
      Map<Integer,OEAtomBase> rgPosToAtomMap = new TreeMap<Integer,OEAtomBase>();
      OEAtomBaseIter atIt = mol.GetAtoms();
      int atMapIdx = 1;
      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         int charge = at.GetFormalCharge();        // charge is RGroup number

         if( charge > 0 && at.GetAtomicNum() == OEElemNo.U )
         {  if( rgPosToAtomMap.containsKey(charge) )
               throw new Error(String.format("U+%d found multiple times", charge));

            // Increment charge so that core smiles contains charge + 100
            charge += 100;
            at.SetFormalCharge(charge);

            rgPosToAtomMap.put(charge,at);

         }else
         {  at.SetMapIdx(atMapIdx++);
            at.SetImplicitHCount(0);
         }
      }
      atIt.delete();

      StringBuilder eductSmi = new StringBuilder();
      // initialize eductSmiles with core smiles
      eductSmi.append(oechem.OECreateSmiString(mol,SMIFlags));

      // loop over rgroups and append disconnected fragments to productSmiles
      for(Entry<Integer, OEAtomBase> posAtom : rgPosToAtomMap.entrySet())
      {  int pos = posAtom.getKey();
         OEAtomBase at = posAtom.getValue();
         at.SetAtomicNum(0);
         at.SetFormalCharge(0);
         at.SetMapIdx(atMapIdx);

         // decrement rgroup charge so that rgroup has correct numbering
         eductSmi.append(".[U+").append(pos-100).append("][*:").append(atMapIdx).append(']');
         atMapIdx++;
      }

      String productSmi = oechem.OECreateSmiString(mol,SMIFlags);
      String smirks =  eductSmi + ">>" + productSmi;
      System.err.println("ScaffoldTransform: " + smirks);

      return smirks;
   }


   public static void main2(String argv[])
   {
      OELibraryGen libgen = new OELibraryGen(
               "[O:1]=[C:2][Cl:3].[N:4][H:5]>>[O:1]=[C:2][N:4]");
      OEGraphMol mol = new OEGraphMol();
      oechem.OEParseSmiles(mol, "CC(=O)Cl");
      libgen.SetStartingMaterial(mol, 0);
      mol.Clear();
      oechem.OEParseSmiles(mol, "NCC");
      libgen.SetStartingMaterial(mol, 1);
      for (OEMolBase product : libgen.GetProducts())
      {
         String smi = oechem.OECreateCanSmiString(product);
         System.out.println("product smiles = " + smi);
      }
   }

   private void delete()
   {  libgen.delete();
   }


   private static void exitWithHelp(String msg, Options options)
   {  System.err.println(msg);
      String myName= "SDFEnumerator";
      HelpFormatter formatter = new HelpFormatter();
      String head = IOUtil.getResource(SDFEnumerator.class, myName+".txt");
      formatter.printHelp( myName, head, options, "", false );

      System.exit(1);
   }
}
