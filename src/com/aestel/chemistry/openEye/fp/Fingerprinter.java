/*
   Copyright 2006-2014 Man-Ling Lee & Alberto Gobbi

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Contact: aestelSW@gmail.com
*/

package com.aestel.chemistry.openEye.fp;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OEMolBase;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;
import openeye.oechem.oemolostream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import com.aestel.utility.IntArrayList;

/**
 * Command line program to compute various kinds of fingerprints and
 * store them as hex encoded string in a tag of the sdf file or to output the
 * fingerprint with and ID in tab separated format.
 *
 * The fingerprinter relies on a mapping file that maps a string representation of
 * the structural feature to the index position in the fingerprint.
 * The couputation of the bit position is done by an implementation of
 * {@link StructureCodeNameIterator}.
 * The mapping onto the bit position is done by the {@link StructureCodeMapper}.
 *
 * This gets about 1M maccs fingerprints into 64MB memory in 1425sec on windows.
 *
 * To recreate the dictionary files for the linear7 linear7*4 fingerprints:
 * - remove existing dictionary file *Map.tab from config/fp.
 * - run with the following options:
 *    -i c:\tmp\allSmdi.ism -o c:\tmp\t.csv -fpType linear7*4 -writeCodeMap
 * - run the FPDictionarySorter with the following options:
 *    -fpType linear7*4 -i c:\tmp\allSmdi.ism -sampleFract .1
 *
 * Note: that after this operation newly created fingerprints will not be
 *       backward compatible.
 *
 * @author albertgo
 *
 */
public class Fingerprinter
{  private final StructureCodeNameIterator sCodeIter;
   private final StructureCodeMapper mapper;


   public Fingerprinter(StructureCodeNameIterator sCodeIter, StructureCodeMapper mapper)
   {  this.sCodeIter = sCodeIter;
      this.mapper = mapper;
   }


   public Fingerprint getFingerprint(OEMolBase mol)
   {  sCodeIter.init(mol);
      IntArrayList bits = new IntArrayList(200);
      while(sCodeIter.hasNext())
      {  String sCodeName = sCodeIter.next();
         int pos = mapper.getIndex(sCodeIter.getType(), sCodeName);
         if(pos >= 0)
            bits.add(pos);
      }
      return new SparseFingerprint(bits.toArray());
   }

   StructureCodeMapper getMapper()
   {  return mapper;
   }

   private void writeDictionary()
   {  if(!(mapper instanceof LearningStrcutureCodeMapper))
         throw new Error("-writeCodeMap may not be used with this mapper!");
      ((LearningStrcutureCodeMapper)mapper).writeDictionary();
   }

   /**
    *
    * @param type type of fingerprints to generate.
    * @param addNewFragments if true new fragments are added to the dictionary use
    *                        this to update or create a new dictionary.
    * @param hashNewFragments hash new fragments into a bit-range exceeding the
    *                         bit-range of the dictionary.
    */
   public static Fingerprinter createFingerprinter(String type,
         boolean addNewFragments, boolean hashNewFragments)
   {  StructureCodeNameIterator generator;
      StructureCodeMapper mapper;
      if ("maccs".equals(type))
      {  generator = SmartsCodeNameIterator.createFromXML(Constants.MACCSSmartsFile);
         if(hashNewFragments) throw new Error("hashNewFragments not supported for maccs");
         mapper = new TABDictionaryStructureCodeMapper("maccsMap.tab", addNewFragments, false);

      } else if ("linear7".equals(type))
      {  generator = new LinearCodeNameGenerator("lin7", 7, 99);
         mapper = new TABDictionaryStructureCodeMapper("linear7Map.tab",
                                             addNewFragments, hashNewFragments);

      } else if ("linear7*4".equals(type))
      {  generator = new LinearCodeNameGenerator("lin74", 7, 4);
         mapper = new TABDictionaryStructureCodeMapper("linear74Map.tab",
                                             addNewFragments, hashNewFragments);

      } else if ("HashLinear7*4".equals(type))
      {  generator = new LinearCodeNameGenerator("HLin74", 7, 4);
         if(hashNewFragments) throw new Error("hashNewFragments not supported for HLin74");
         if(addNewFragments)  throw new Error("addNewFragments not supported for HLin74");
         mapper = new HashStructureCodeMapper(0, 16348);

      } else
      {  throw new Error("Unknown fingerprint type: " + type);
      }

      return new Fingerprinter(generator, mapper);
   }


   public void close()
   {  if(mapper != null)    mapper.close();
      if(sCodeIter != null) sCodeIter.close();
   }


   public static void main(String...args) throws IOException
   {  long start = System.currentTimeMillis();
      long iCounter = 0;

      // create command line Options object
      Options options = new Options();
      Option opt = new Option("in",true, "input file [.ism,.sdf,...]");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("out",true, "output file .tsv or oe-supported");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("idTag",true, "field with ID (default title)");
      opt.setRequired(false);
      options.addOption(opt);

      opt = new Option("fpType",true,
         "fingerPrintType: maccs|linear7|linear7*4|HashLinear7*4\n"
        +"   maccs: generate maccs keys\n"
        +"   linear7 generate 7 bonds long linear fingerprints (210k known rest hashed)\n"
        +"   linear7*4 linear 7 bonds if more than 4 atoms code atoms as * (5.1k known rest hashed)\n"
        +"   HashLinear7*4: as linear7*4 but hashed to 16k");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("format",true,
          "folded512|folded2048|bitList|fragList|none\n"
         +"   folded512/2048: hex encoded 512/2048 bits\n"
         +"   bitList: list of bitpositions\n"
         +"   none: no fp output for use with writeCodeMap");
      opt.setRequired(true);
      options.addOption(opt);

      opt = new Option("writeCodeMap",false,
                        "Overwrite the codeMap file at the end of processing");
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

      String idTag = null;
      if(cmd.hasOption("idTag"))
         idTag = cmd.getOptionValue("idTag");

      String outformat = cmd.getOptionValue("format").toLowerCase().intern();
      if(args.length != 0)
      {  exitWithHelp(options);
      }

      String type = cmd.getOptionValue("fpType");
      boolean updateDictionaryFile = cmd.hasOption("writeCodeMap");
      boolean hashUnknownFrag = true;
      if( type.equals("HashLinear7*4") ) hashUnknownFrag = false;
      if( type.equals("maccs") )         hashUnknownFrag = false;
      if( updateDictionaryFile )         hashUnknownFrag = false;
      Fingerprinter fprinter = createFingerprinter(type,
                                          updateDictionaryFile, hashUnknownFrag);
      OEMolBase mol = new OEGraphMol();

      String inFile  = cmd.getOptionValue("in");
      String outFile = cmd.getOptionValue("out");
      oemolistream ifs = new oemolistream(inFile);

      Runtime rt = Runtime.getRuntime();
      Outputter out;
      if( outFile.endsWith(".txt") || outFile.endsWith(".tab"))
         out = new TabOutputter(fprinter.getMapper(), outFile, outformat);
      else
         out = new OEOutputter(fprinter.getMapper(), outFile, type, outformat);

      while(oechem.OEReadMolecule(ifs, mol))
      {  iCounter++;
         Fingerprint fp = fprinter.getFingerprint(mol);

         String id;
         if(idTag == null)
            id = mol.GetTitle();
         else
            id = oechem.OEGetSDData(mol, idTag);

         if(iCounter % 100 == 0) System.err.print(".");
         if(iCounter % 4000 == 0)
         {  System.err.printf( " %d %dsec\tt=%d f=%d u=%d m=%d tf=%d\n",
                  iCounter, (System.currentTimeMillis()-start)/1000,
                  rt.totalMemory()/1024, rt.freeMemory()/1024,
                  (rt.totalMemory() - rt.freeMemory())/1024, rt.maxMemory()/1024,
                  (rt.freeMemory() + (rt.maxMemory() - rt.totalMemory())) / 1024 );
         }

         out.output(id, mol, fp);
      }

      System.err.printf("Fingerprinter: Read %d structures in %d sec\n",
            iCounter, (System.currentTimeMillis()-start)/1000);

      if(updateDictionaryFile ) fprinter.writeDictionary();
      out.close();
      fprinter.close();
   }

   private static void exitWithHelp(Options options) {
      HelpFormatter formatter = new HelpFormatter();
      formatter.printHelp( "fingerprinter", options );
      System.exit(1);
   }
}

interface Outputter
{  void output(String id, OEMolBase mol, Fingerprint pf);
   void close();
}

class OEOutputter implements Outputter
{  oemolostream ofs;
   private String outFormat;
   private final StringBuilder sb = new StringBuilder(500);
   private String tagName;
   private final StructureCodeMapper mapper;

   OEOutputter(StructureCodeMapper mapper, String fName, String fpType, String outFormat)
      throws FileNotFoundException
   {  ofs = new oemolostream(fName);
      this.outFormat = outFormat;
      this.tagName = fpType + "_" + outFormat;
      this.mapper = mapper;
   }

   @Override
   public void output(String id, OEMolBase mol, Fingerprint fp)
   {  if(outFormat.equalsIgnoreCase("folded512"))
         oechem.OESetSDData(mol, tagName, fp.fold(512).getHexString());
      else if(outFormat.equalsIgnoreCase("folded2048"))
         oechem.OESetSDData(mol, tagName, fp.fold(2048).getHexString());
      else if(outFormat.equalsIgnoreCase("bitlist"))
      {  sb.setLength(0);
         sb.append(id).append('\t');
         for(int idx :fp.getBits())
            sb.append(idx).append(',');
         sb.setLength(sb.length()-1);

         oechem.OESetSDData(mol, tagName, sb.toString());
      }
      else if(outFormat.equalsIgnoreCase("fragList"))
      {  sb.setLength(0);
         sb.append(id).append('\t');
         for(int idx :fp.getBits())
            sb.append(mapper.getName(idx)).append(',');
         sb.setLength(sb.length()-1);

         oechem.OESetSDData(mol, tagName, sb.toString());
      }

      oechem.OEWriteMolecule(ofs, mol);
   }

   @Override
   public void close()
   {  ofs.close();
   }
}

class TabOutputter implements Outputter
{  PrintStream out;
   private String outFormat;
   private final StringBuilder sb = new StringBuilder(500);
   private final StructureCodeMapper mapper;

   TabOutputter(StructureCodeMapper mapper, String fName, String outFormat) throws FileNotFoundException
   {  PrintStream tmpout = System.out;
      if(! (".tab".equalsIgnoreCase(fName) || ".txt".equalsIgnoreCase(fName)) )
         tmpout = new PrintStream(new FileOutputStream(fName));
      out = tmpout;
      this.outFormat = outFormat;
      this.mapper = mapper;
   }

   @Override
   public void output(String id, OEMolBase mol, Fingerprint fp)
   {  if(outFormat.equalsIgnoreCase("folded512"))
         out.println(id + "\t" + fp.getNBits()+ "\t" + fp.fold(512).getHexString());
      else if(outFormat.equalsIgnoreCase("folded2048"))
         out.println(id + "\t" + fp.getNBits()+ "\t" + fp.fold(2048).getHexString());
      else if(outFormat.equalsIgnoreCase("bitlist"))
      {  sb.setLength(0);
         sb.append(id).append('\t');
         for(int idx :fp.getBits())
            sb.append(idx).append(',');
         sb.setLength(sb.length()-1);

         out.println(sb);
      } else if(outFormat.equalsIgnoreCase("fragList"))
      {  sb.setLength(0);
         sb.append(id).append('\t');
         for(int idx :fp.getBits())
            sb.append(mapper.getName(idx)).append(',');
         sb.setLength(sb.length()-1);

         out.println(sb);
      }
   }

   @Override
   public void close()
   {  out.close();
   }
}
