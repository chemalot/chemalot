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

package com.genentech.application.property;

/*
Calculation of polar surface area based on fragment contributions (TPSA)
Peter Ertl, Novartis Pharma AG, Basel, August 2000
peter.ertl@pharma.novartis.com || peter.ertl@altavista.net
The program uses the SMILES Daylight Toolkit module.

For description of the methodology see :
P. Ertl, B. Rohde, P. Selzer
Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
Contributions and Its Application to the Prediction of Drug Transport
Properties, J.Med.Chem. 43, 3714-3717, 2000

The program reads a list of SMILES strings from stdin, and writes SMILES
and the calculated TPSA on stdout.
The program expects "problematic" nitrogens in pentavalent form, i.e. nitro
groups as -N(=O)=O, N-oxides as c1ccccn1=O, and azides as -N=N#N.

This contrib program comes with absolutely no warranty of any kind.

Rev: 29 Mar 2001

 */


import openeye.oechem.*;
import java.util.*;
import org.apache.commons.cli.*;

public class TPSA {

   public double calculateTPSA(OEGraphMol mol, boolean countP, boolean countS) {

      double psa = 0.0;
      List<String> errArray = new ArrayList<String>(); // array of error messages

      // loop through all atoms
      for (OEAtomBaseIter aiter = mol.GetAtoms(); aiter.hasNext();) {
         OEAtomBase atom = aiter.next();

         // Check atom type
         int an = atom.GetAtomicNum();

         // N, O, P and S are consider polar, can treat S and P as non-polar by using the
         //	countP and countS flags
         if (an != OEElemNo.N && an != OEElemNo.O && an != OEElemNo.P
               && an != OEElemNo.S )
            continue;

         int nh = atom.GetTotalHCount();
         int q = atom.GetFormalCharge();

         // Check if atom is in 3 member ring
         boolean isIn3ring = oechem.OEAtomIsInRingSize(atom, 3);

         // Find the number of -, +, #, : bonds originating from the atom
         int nsingle = 0;
         int ndouble = 0;
         int ntriple = 0;
         int naromatic = 0;

         // Number of neighbors
         int nv = 0;
         for (OEBondBaseIter biter = atom.GetBonds(); biter.hasNext();) {
            // Count the number of each type of bond
            OEBondBase bond = biter.next();
            /* do not count bonds that are bonded to H */
            if ( (bond.GetBgn().GetAtomicNum() == OEElemNo.H) || (bond.GetEnd().GetAtomicNum() == OEElemNo.H) ){
               continue;
            }
            nv++;

            int bondType = bond.GetOrder();
            //System.out.println("bondType " + bondType);
            if (bond.IsAromatic() ) {
               //System.out.println("aromatic bond " + bond.GetBgnIdx() + "-" + bond.GetEndIdx());
               naromatic++;
            } else {
               switch (bondType) {
               case 1:
                  nsingle++;
                  break;
               case 2:
                  ndouble++;
                  break;
               case 3:
                  ntriple++;
                  break;
               }
            }
         } // end for (OEBondBaseIter

         // Find the PSA contribution for this fragment (atom with hydrogens)
         double p = -1.0;

         switch (an) {
         case OEElemNo.N:
            if (nv == 1) {
               if (nh == 0 && ntriple == 1 && q == 0) // N#
                  p = 23.79;
               else if (nh == 1 && ndouble == 1 && q == 0) // [NH]=
                  p = 23.85;
               else if (nh == 2 && nsingle == 1 && q == 0) // [NH2]-
                  p = 26.02;
               else if (nh == 2 && ndouble == 1 && q == 1) // [NH2+]=
                  p = 25.59;
               else if (nh == 3 && nsingle == 1 && q == 1) // [NH3+]-
                  p = 27.64;
               else
                  errArray.add("Standard TPSA parameter not found for atom N" + atom.GetIdx());
            } else if (nv == 2) {
               if (nh == 0 && nsingle == 1 && ndouble == 1 && q == 0) // =N-
                  p = 12.36;
               else if (nh == 0 && ntriple == 1 && ndouble == 1 && q == 0) // =N#
                  p = 13.60;
               else if (nh == 1 && nsingle == 2 && q == 0 && !isIn3ring) // -[NH]-
                  p = 12.03;
               else if (nh == 1 && nsingle == 2 && q == 0 && isIn3ring) // -[NH]-r3
                  p = 21.94;
               else if (nh == 0 && ntriple == 1 && nsingle == 1 && q == 1) // -[N+]#
                  p = 4.36;
               else if (nh == 1 && ndouble == 1 && nsingle == 1 && q == 1) // -[NH+]=
                  p = 13.97;
               else if (nh == 2 && nsingle == 2 && q == 1) // -[NH2+]-
                  p = 16.61;
               else if (nh == 0 && naromatic == 2 && q == 0) // :[n]:
                  p = 12.89;
               else if (nh == 1 && naromatic == 2 && q == 0) // :[nH]:
                  p = 15.79;
               else if (nh == 1 && naromatic == 2 && q == 1) // :[nH+]:
                  p = 14.14;
               else
                  errArray.add("Standard TPSA parameter not found for atom N" + atom.GetIdx());

            } else if (nv == 3) {
               if (nh == 0 && nsingle == 3 && q == 0 && !isIn3ring) // -N(-)-
                  p = 3.24;
               else if (nh == 0 && nsingle == 3 && q == 0 && isIn3ring) // -N(-)-r3
                  p = 3.01;
               else if (nh == 0 && nsingle == 1 && ndouble == 2 && q == 0) // -N(=)=
                  p = 11.68;
               else if (nh == 0 && nsingle == 2 && ndouble == 1 && q == 1) // =[N+](-)-
                  p = 3.01;
               else if (nh == 1 && nsingle == 3 && q == 1) // -[NH+](-)-
                  p = 4.44;
               else if (nh == 0 && naromatic == 3 && q == 0) // :[n](:):
                  p = 4.41;
               else if (nh == 0 && nsingle == 1 && naromatic == 2 && q == 0) // -:[n](:):
                  p = 4.93;
               else if (nh == 0 && ndouble == 1 && naromatic == 2 && q == 0) // =:[n](:):
                  p = 8.39;
               else if (nh == 0 && naromatic == 3 && q == 1) // :[n+](:):
                  p = 4.10;
               else if (nh == 0 && nsingle == 1 && naromatic == 2 && q == 1) // -:[n+](:):
                  p = 3.88;
               else
                  errArray.add("Standard TPSA parameter not found for atom N" + atom.GetIdx());
            } else if (nv == 4) {
               if (nh == 0 && nsingle == 4 && q == 1) // -[N+](-)(-)-
                  p = 0.00;
               else
                  errArray.add("Standard TPSA parameter not found for atom N" + atom.GetIdx());
            }
            if (p < 0.) { // N with non-standard valency
               errArray.add("Non-standard valency for atom N" + atom.GetIdx());
               p = 30.5 - nv * 8.2 + nh * 1.5;
               if (p < 0.)
                  p = 0.;
            }
            break;

         case OEElemNo.O:
            if (nv == 1) {
               if (nh == 0 && ndouble == 1 && q == 0) // O=
                  p = 17.07;
               else if (nh == 1 && nsingle == 1 && q == 0) // [OH]-
                  p = 20.23;
               else if (nh == 0 && nsingle == 1 && q == -1) // [O-]-
                  p = 23.06;
               else
                  errArray.add("Standard TPSA parameter not found for atom O" + atom.GetIdx());
            } else if (nv == 2) {
               if (nh == 0 && nsingle == 2 && q == 0 && !isIn3ring) // -O-
                  p = 9.23;
               else if (nh == 0 && nsingle == 2 && q == 0 && isIn3ring) // -O-r3
                  p = 12.53;
               else if (nh == 0 && naromatic == 2 && q == 0) // :o:
                  p = 13.14;
               else
                  errArray.add("Standard TPSA parameter not found for atom O" + atom.GetIdx());
            }
            if (p < 0.) { // O with non-standard valency
               p = 28.5 - nv * 8.6 + nh * 1.5;
               if (p < 0.)
                  p = 0.;
            }

            break;

         case OEElemNo.P:
            if (countP){
               if (naromatic > 0)
                  p = 0.0;
               else if (nv == 2 && nsingle == 1 && ndouble == 1) // [P](-*)=*
                  p = 34.14;
               else if (nv == 3 && nsingle == 2 && ndouble == 1) // [P](=*)(-*)-*
                  p = 23.47;
               else if (nv == 3 && nsingle == 3) // [P](-*)(-*)-*
                  p = 13.59;
               else if (nv == 4 && nsingle == 3 && ndouble == 1) // [P](-*)(-*)(-*)=*
                  p = 9.81;
               else
                  errArray.add("Standard TPSA parameter not found for atom P" + atom.GetIdx());
            } else{
               p=0.0;
            }

            break;

         case OEElemNo.S:
            if (countS) {
               if (nv == 1 ) {
                  if (nh == 0 && ndouble == 1 && q == 0) // S=
                     p = 32.09;
                  else if (nh == 1 && nsingle == 1 && q == 0) // [SH]-
                     p = 38.80;
                  else if (nh == 0 && nsingle == 1 && q == -1) // [S-]-
                     p = 38.80;
                  else
                     errArray.add("Standard TPSA parameter not found for atom S" + atom.GetIdx());
               }
               else if (nv == 2 && naromatic == 2) // [s](:*):*
                  p = 28.24;
               else if (nv == 2 && nsingle == 2) // [S](-*)-*
                  p = 25.30;
               else if (nv == 3 && naromatic == 2 && ndouble == 1) // [s](:*)(:*)=*
                  p = 21.70;
               else if (nv == 3 && nsingle == 2 && ndouble == 1) // [S](-*)(-*)=*
                  p = 19.21;
               else if (nv == 4 && nsingle == 2 && ndouble == 2) // [S](-*)(-*)(=*)=*
                  p = 8.38;
               else
                  errArray.add("Standard TPSA parameter not found for atom S" + atom.GetIdx());
            }else {
               p=0.0;
            }

         }
         psa += p;
         //System.out.print(oechem.OEGetAtomicSymbol(an) + atom.GetIdx() + " " + p + ", ");
      }

      Iterator<String> iter = errArray.iterator();
      while (iter.hasNext()) {
         System.err.println(mol.GetTitle() + " " + iter.next());
      }
      return psa;
   }

   public static void main(String args[]) {
      try {
         boolean countP = false;
         boolean countS = false;
         String[] molFiles;
         // create Options object
         Options options = new Options();

         // add  options
         options.addOption("P", false, "Count phosphorus atoms");
         options.addOption("S", false, "Count sulfur atoms");

         CommandLineParser parser = new PosixParser();
         CommandLine cmd = parser.parse(options, args);
         if (cmd.hasOption("P"))	countP = true;
         if (cmd.hasOption("S")) countS = true;

         //get rest of arguments as an array of string
         molFiles = cmd.getArgs();
         if (molFiles.length < 1) {
            // automatically generate the help statement
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "TPSA [options] <mol files> ... ", options );
            System.exit(1);
         }

         TPSA myTPSA = new TPSA();
         myTPSA.runTest(countP, countS);


         OEGraphMol mol = new OEGraphMol();
         oemolistream ifs = new oemolistream();;


         for( int i=0; i< molFiles.length; i++) {
            ifs= new oemolistream(molFiles[i]);
            while (oechem.OEReadMolecule(ifs, mol)) {

               double tpsa = myTPSA.calculateTPSA(mol, countP, countS);
               if (tpsa < 0) {
                  System.out.printf("%s 0\n", mol.GetTitle());
               } else {
                  System.out.printf("%s %d\n", mol.GetTitle(), (int)tpsa);
               }

            }
         }
         ifs.close();
      } catch (ParseException e){
         e.printStackTrace();
      } catch (Exception e) {
         e.printStackTrace();
      }
   }


   /* unit test code */
   private void runTest( boolean countP, boolean countS) {
      String canSmi;
      OEGraphMol mol = new OEGraphMol();
      LinkedHashMap<String,String> smilesMap = new LinkedHashMap<String,String>();
      loadHashMap(smilesMap);
      Iterator<String> iter = smilesMap.keySet().iterator();
      oemolostream ofs = new oemolostream();
      ofs.open("TPSA_test.sdf");
      while (iter.hasNext()) {
         String key = (String) iter.next();
         String PSA = (String) smilesMap.get(key);
         mol.Clear();
         if (oechem.OEParseSmiles(mol, key)) {
            mol.SetTitle(key);
            System.out.print(key + ": ");
            double tpsa = calculateTPSA(mol, countP, countS);
            mol.SetTitle(key);
            oechem.OEAddSDData(mol, "TPSA", String.valueOf(tpsa));
            canSmi = oechem.OECreateCanSmiString(mol);
            oechem.OEAddSDData(mol, "CanSmiles", canSmi);
            System.out.printf("Lit: %s\tTPSA: %.2f\n", PSA, tpsa);
            oechem.OEWriteMolecule(ofs, mol);
         }else {
            System.err.println("Error parsing the SMILES string: " + key);
         }
      }
   }

   /* basic test cases */
   private void loadHashMap(Map <String,String> smilesMap){

      smilesMap.put("CN(C)C", "3.24");
      smilesMap.put("CN=C", "12.36");
      smilesMap.put("C#N", "23.79");
      smilesMap.put("CN(=O)=O", "11.68");
      smilesMap.put("C=N#C", "13.60");
      smilesMap.put("CN1CC1", "3.01");
      smilesMap.put("CNC", "12.03");
      smilesMap.put("C1CN1", "21.94");
      smilesMap.put("C=N", "23.85");
      smilesMap.put("CN", "26.02");
      smilesMap.put("C[N+](C)(C)C", "0.00");
      smilesMap.put("C[N+](=C)C", "3.01");
      smilesMap.put("C[N+]#[C-]", "4.36");
      smilesMap.put("C[NH+](C)C", "4.44");
      smilesMap.put("C[NH+]=C", "13.97");
      smilesMap.put("C[NH2+]C", "16.61");
      smilesMap.put("C=[NH2+]", "25.59");
      smilesMap.put("C[NH3+]", "27.64");
      smilesMap.put("c1ccncc1", "12.89");
      smilesMap.put("c1ccn2cccc2c1", "4.41");
      smilesMap.put("Cn1cccc1", "4.93");
      smilesMap.put("c1ccn(=C)cc1", "8.39");
      smilesMap.put("c1cc[nH]c1", "15.79");
      smilesMap.put("c1cc[n+]2ccccc2c1", "4.10");
      smilesMap.put("C[n+]1ccccc1", "3.88");
      smilesMap.put("c1cc[nH+]cc1", "14.14");
      smilesMap.put("COC", "9.23");
      smilesMap.put("C1CO1", "12.53");
      smilesMap.put("C=O", "17.07");
      smilesMap.put("OC", "20.23");
      smilesMap.put("C[O-]", "23.06");
      smilesMap.put("c1ccoc1", "13.14");
      smilesMap.put("CSC", "25.30");
      smilesMap.put("C=S", "32.09");
      smilesMap.put("CS(=C)C", "19.21");
      smilesMap.put("CS(=O)(=O)C", "8.38");
      smilesMap.put("CS", "38.80");
      smilesMap.put("c1ccsc1", "28.24");
      smilesMap.put("C=s1cccc1", "21.70");
      smilesMap.put("CP(C)C", "13.59");
      smilesMap.put("CP=C", "34.14");
      smilesMap.put("OP(=O)(O)O", "9.81");
      smilesMap.put("CP(=O)O", "23.47");

      /*selected structures from table 3 / fig 3 of Ertl et al. */
      smilesMap.put("OC(P(O)(O)=O)=O", "94.8"); /*foscarnet */

      /* SMILES in tabel 4 of Ertl et al. */
      smilesMap.put("O=C([O-])c1ccccc1O", "60.36");
      smilesMap.put("CC(=O)Oc1ccccc1C(=O)[O-]", "66.43");
      smilesMap.put("OCC(O)C(O)C(O)C(O)CO", "121.37");
      smilesMap.put("CC(C)(C)NCC(O)c1cc(O)cc(O)c1", "72.71");
      smilesMap.put("C=CCc1ccccc1OCC(O)CNC(C)C", "41.49");
      smilesMap.put("CC(C)NCC(O)COc1cccc2ccccc12", "41.49");
      smilesMap.put("CC(=O)Nc1ccc(OCC(O)CNC(C)C)cc1", "70.58");
      smilesMap.put("CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", "84.58");
      smilesMap.put("COCCc1ccc(OCC(O)CNC(C)C)cc1", "50.72");
      smilesMap.put("CC12CCC(=O)C=C1CCC3C2CCC4(C)C(O)CCC34", "37.30");
      smilesMap.put("O=C(O)c2cc(N=Nc1ccc(O)c(C(=O)O)c1)ccc2O", "139.78");
      smilesMap.put("CC(=O)CC(c1ccccc1)c3c(O)c2ccccc2oc3=O", "67.51");
      smilesMap.put("CC14CCC(=O)C=C1ccc3c2ccc(C(=O)CO)C2(C)CC(O)C34", "74.60");
      smilesMap.put("CC14CCC(=O)C=C1ccc3c2ccc(O)(C(=O)CO)C2(C)CC(O)C34", "94.83");
      smilesMap.put("CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c2cccc(C1)c2C1", "64.64");
      smilesMap.put("CC4CC3C2CCC1=CC(=O)C=CC1(C)C2(F)C(O)CC3(C)C4(O)C(=O)CO", "94.83");
      smilesMap.put("O=C(O)c3cc(N=Nc2ccc(S(=O)(=O)Nc1ccccn1)cc2)ccc3O", "141.32");


      /*selected SMILES from tabel 5 of Ertl et al */
      smilesMap.put("CNC(=CN(=O)=O)NCCSCc1ccc(CN(C)C)o1", "86.26");
      smilesMap.put("CNC(=CC#N)NCCSCc1csc(N=C(N)N)n1", "125.15");
      smilesMap.put("NC(N)=Nc2nc(c1ccccc1)cs2", "77.3");


   }

}
