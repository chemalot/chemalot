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
package com.genentech.chemistry.openEye.cats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import openeye.oechem.*;

import com.aestel.chemistry.molecule.Atom;
import com.aestel.utility.DataFormat;

public class CATSIndexer
{
//   static final String BASICGroup  = "[$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))]),$([N,n;X2;+0])]";
//   static final String ACIDICGroup = "[$([C,S](=[O,S,P])-[O;H1])]";
//   static final String HDonor      = "[$([N;!H0;v3]),$([N;!H0;+1;v4]),$([O,S;H1;+0]),$([n;H1;+0])] ";
//   static final String HAcceptor   = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]";
//
   static final String BASICGroup  = "[N&!$(N[a])&!$(N*=O)&!$(N~[O,N])&!$(N#*)]";
   static final String ACIDICGroup = "[$([OH]*=O),$(O=*[OH]),$([nH](n)nn),$(n(n)[nH]n),$([NH](C=O)*=O),#8-]";
   static final String HDonor      = "[#7&!H0&!$(n(n)nn),OH&!$(O*=O),SH]";
   static final String HAcceptor   = "[#7X2,$(N#*),#8]";
   static final String CARBONYL    = "C=O";
   static final String URANIUM     = "[U]";

   static final OESubSearch CARBONYLSS = new OESubSearch(CARBONYL);

   String[] allDescriptors = getAllDescriptorNames();

   private final String tagPrefix;
   private final AtomTyperInterface[] myTypes;
   private final int maxBondDist;


   public static final AtomTyperInterface[] typers =
      {  new SSTyper(HDonor,0, "D"),
         new SSTyper(HAcceptor, 1, "A"),
         new SSTyper(BASICGroup, 2, "P"),
         new SSTyper(ACIDICGroup, 3, "N"),
         new LipophilicTyper(4)
      };

   public static final AtomTyperInterface[] rgroupTypers =
      {  typers[0], typers[1], typers[2], typers[3], typers[4],
         new SSTyper(URANIUM, 5, "U")
      };

   // Method to load a tab delimited SMARTS<tab>ID\n list of features for indexing
   public static AtomTyperInterface[] featuresFromFile(String featureFilename) throws IOException
   {
      ArrayList<AtomTyperInterface> featureList = new ArrayList<AtomTyperInterface>();
      if( featureFilename == null || featureFilename.length() == 0 )
         throw new FileNotFoundException("Invalid feature File Name: null");

      BufferedReader featReader = new BufferedReader( new FileReader( featureFilename ) );
      String line = "";

      int featCount = 0;
      while( ( line = featReader.readLine() ) != null )
      {
         if( line.length() < 1 || line.startsWith( "#" ) || line.startsWith( "\"#" ) )
            continue;

         String[] items = line.split( "\t" ); //SMARTS<tab>FEAT_NAME

         String featureSMARTS = items[0].trim();
         String featureName   = items[1].trim();

         featureList.add(new SSTyper(featureSMARTS, featCount, featureName));
         featCount++;
      }

      featReader.close();
      System.err.printf( "Read %d user features\n", featCount );

      return featureList.toArray(new AtomTyperInterface[featureList.size()]);
   }


   private final OEGraphMol mCopy = new OEGraphMol();

   public static enum Normalization
   {  Counts,
      CountsPerAtom,
      CountsPerFeature
   }

   public CATSIndexer(AtomTyperInterface[] myTypes, String tagPrefix, int maxBondDistance)
   {  this.myTypes = myTypes;
      this.tagPrefix = tagPrefix + "CATS";
      this.maxBondDist = maxBondDistance;
   }


   public void compute2DCats(OEMolBase mol, EnumSet<Normalization> normMeth)
   {  mCopy.Clear();
      oechem.OEAddMols(mCopy, mol);

      oechem.OEAssignAromaticFlags(mCopy);
      oechem.OEPerceiveChiral(mCopy);
      oechem.OEAssignHybridization(mCopy);
      oechem.OESuppressHydrogens(mCopy);

      List<OEAtomBase> typedAtoms = new ArrayList<OEAtomBase>();
      List<List<AtomTyperInterface>> atTypeList = new ArrayList<List<AtomTyperInterface>>();
      computeAtomTypes(mCopy, typedAtoms, atTypeList);

      // compute distance and pair count
      int[][][] pairCount = new int[myTypes.length][myTypes.length][maxBondDist+1];
      int[] typeCount = new int[myTypes.length];
      compute2DPairCount(typedAtoms, atTypeList, pairCount, typeCount);

      //System.err.println(distArrayToString(pairCount));
      // normalize as in:
      // Fechner U, Schneider G: QSAR Comb Sci 2004, 23:19 22.

      if( normMeth.contains(Normalization.Counts) )
      {  String myTagPrefix = tagPrefix + "2C_";
         for(int i=0; i<myTypes.length; i++)
            for(int j=0; j<=i; j++)
               for(int d=0; d<=maxBondDist; d++)
                  oechem.OESetSDData(mol,
                      myTagPrefix + d + "_" + myTypes[i].getTypeName() + myTypes[j].getTypeName(),
                      Integer.toString(pairCount[i][j][d]));
      }

      if( normMeth.contains(Normalization.CountsPerAtom) )
      {  String myTagPrefix = tagPrefix + "2A_";
         double atCount = getAtCount(mol);
         for(int i=0; i<myTypes.length; i++)
            for(int j=0; j<=i; j++)
               for(int d=0; d<=maxBondDist; d++)
                  oechem.OESetSDData(mol,
                      myTagPrefix + d + "_" + myTypes[i].getTypeName() + myTypes[j].getTypeName(),
                      DataFormat.round(pairCount[i][j][d]/atCount, 3));
      }

      if( normMeth.contains(Normalization.CountsPerFeature) )
      {  String myTagPrefix = tagPrefix + "2F_";
         for(int i=0; i<myTypes.length; i++)
            for(int j=0; j<=i; j++)
               for(int d=0; d<=maxBondDist; d++)
               {  double val = pairCount[i][j][d];
                  if ( val > 0 )
                     val = val / (typeCount[i] + typeCount[j]);

                  oechem.OESetSDData(mol,
                      myTagPrefix + d + "_" + myTypes[i].getTypeName() + myTypes[j].getTypeName(),
                      DataFormat.round(val,3));
               }
      }
   }

   public void printDistanceDescriptors(PrintStream out, OEMolBase mol)
   {  mCopy.Clear();
      oechem.OEAddMols(mCopy, mol);

      oechem.OEDetermineConnectivity(mCopy);
      oechem.OEFindRingAtomsAndBonds(mCopy);
      oechem.OEPerceiveBondOrders(mCopy);
      oechem.OEAssignImplicitHydrogens(mCopy);
      oechem.OEAssignFormalCharges(mCopy);

      oechem.OEAssignAromaticFlags(mCopy);
      oechem.OEPerceiveChiral(mCopy);
      oechem.OEAssignHybridization(mCopy);
      oechem.OESuppressHydrogens(mCopy);

      OEAtomBase[] atoms = getAtoms(mCopy.GetAtoms());
      String smi = oechem.OECreateSmiString(mCopy);

      Map<Path,String> namedPaths = new HashMap<Path, String>();
      OEMatchBaseIter mIt = CARBONYLSS.Match(mCopy);
      while(mIt.hasNext())
      {  OEMatchBase m = mIt.next();
         Path p = new Path(m);
         namedPaths.put(p, "transAmide");
      }
      mIt.delete();

      for(int i=0; i<atoms.length; i++)
      {  for(int j=0; j<i; j++)
         {  PathDescription desc = new PathDescription(mCopy, atoms[i], atoms[j]);
            if( desc.getPathLength() > maxBondDist ) break;
            desc.computeDescription();

            out.print(smi);
            out.print('\t');
            out.print(desc.getSmiles());
            out.print('\t');
            out.print(DataFormat.round(oechem.OEGetDistance(mCopy, atoms[i], atoms[j]), 2));

            for(String name: allDescriptors)
            {  out.print('\t');
               out.print(DataFormat.round(desc.getDescriptor(name),3));
            }

            out.println();
            desc.close();
         }
      }
   }


   public void printDescriptorHeader(PrintStream out)
   {  out.print("SMI\tPathSMI\tdist3d");
      for(String name: allDescriptors)
      {  out.print('\t');
         out.print(name);
      }
      out.println();
   }

   /**
    * This also deletes the atIt!
    */
   static OEAtomBase[] getAtoms(OEAtomBaseIter atIt)
   {  List<OEAtomBase> atomL = new ArrayList<OEAtomBase>();

      while(atIt.hasNext())
      {  OEAtomBase at = atIt.next();
         atomL.add(at);
      }
      atIt.delete();

      return atomL.toArray(new OEAtomBase[atomL.size()]);
   }


   /**
    * @param typedAtoms list of atoms with types
    * @param atTypeList list of types for each atom
    * @param pairCount filled with count of ouccurence index by atom2Idx,atom2Idx,distanceBin
    * @param typeCount filled with occerence count of each atomType ove whole molecule
    */
   private void compute2DPairCount(List<OEAtomBase> typedAtoms, List<List<AtomTyperInterface>> atTypeList,
            int[][][] pairCount,
            int[] typeCount)
   {
      for(int i=0; i<typedAtoms.size(); i++)
      {  OEAtomBase at1 = typedAtoms.get(i);
         List<AtomTyperInterface> types1 = atTypeList.get(i);

         for( AtomTyperInterface type : types1)
            typeCount[type.getTypeIdx()]++;

         for(int j=0; j<=i; j++ )
         {  OEAtomBase at2 = typedAtoms.get(j);
            List<AtomTyperInterface> types2 = atTypeList.get(j);

            int dist;
            if( i == j )
               dist = 0;
            else
            {  dist = oechem.OEGetPathLength(at1,at2, maxBondDist);
               if( dist == 0 ) continue;
            }

            // record paircount
            for( AtomTyperInterface t1: types1 )
            {  for( AtomTyperInterface t2: types2 )
               {  // order types
                  if( t1.getTypeIdx() < t2.getTypeIdx() )
                     pairCount[t2.getTypeIdx()][t1.getTypeIdx()][dist]++;
                  else
                     pairCount[t1.getTypeIdx()][t2.getTypeIdx()][dist]++;
               }
            }
         }
      }
   }


   /**
    * @param typedAtoms filled with list of atoms that have atom types
    * @param atTypeList filled with list of AtomTyperInterface of each atom
    */
   private void computeAtomTypes(OEMolBase mol,
                                 List<OEAtomBase> typedAtoms,
                                 List<List<AtomTyperInterface>> atTypeList)
   {
      OEAtomBaseIter atit = mol.GetAtoms();

      // get list of typed atoms and their types
      while(atit.hasNext())
      {  OEAtomBase at = atit.next();
         ArrayList<AtomTyperInterface> atTypes = new ArrayList<AtomTyperInterface>();
         for(AtomTyperInterface type:myTypes)
            if( type.isType(at) )   atTypes.add(type);
         if( atTypes.size() > 0 )
         {  typedAtoms.add(at);
            atTypeList.add(atTypes);

//            System.err.printf("%s %s\n", com.genentech.oechem.tools.Atom.getAtomName(at),
//                                         toString(atTypes));
         }
      }
      atit.delete();
   }


   private static int getAtCount(OEMolBase mol)
   {  OEAtomBaseIter atit = mol.GetAtoms();
      int cnt = 0;
      while(atit.hasNext())
      {  OEAtomBase at = atit.next();
         if( at.GetAtomicNum() > 1 ) cnt ++;
      }
      atit.delete();
      return cnt;
   }


   public void close()
   {  for(AtomTyperInterface type: myTypes)
         type.close();
      mCopy.delete();
   }

   public static String toString(ArrayList<AtomTyperInterface> types)
   {  StringBuilder sb = new StringBuilder();
      for(AtomTyperInterface t: types)
         sb.append(t.getTypeName()).append(", ");

      if( sb.length() > 2 ) return sb.substring(0, sb.length()-2);
      return "";
   }

   public String distArrayToString(int[][][] pairCount)
   {  StringBuilder sb = new StringBuilder();
      for(int i=0; i<pairCount.length; i++)
      {  sb.append(myTypes[i].getTypeName()).append(":\t");
         for(int j = 0; j<pairCount[i].length; j++)
         {  sb.append(myTypes[j].getTypeName()).append(':');
            for( int k = 0; k<pairCount[i][j].length; k++)
               sb.append(pairCount[i][j][k]).append('+');
            if('+' == sb.charAt(sb.length()-1)) sb.setLength(sb.length()-1);
            sb.append("\t");
         }
         sb.append("\n");
      }
      return sb.toString();
   }


   private String[] getAllDescriptorNames()
   {  List<String> dNames = new ArrayList<String>();
      dNames.add("pathLen"      );
      dNames.add("vdwSum"       );
      dNames.add("bondOrderSum" );
      dNames.add("totVal"       );
      dNames.add("totSubst"     );
      dNames.add("nRingBonds"   );
      dNames.add("nSingle"      );
      dNames.add("nDouble"      );
      dNames.add("nCisBond"     );
      dNames.add("nTransBond"   );
      dNames.add("nTripple"     );
      dNames.add("nArom"        );

      for(int i=1; i<=maxBondDist; i++)
      {  dNames.add("atVal"+i);
         dNames.add("atSubst"+i);
         dNames.add("atIsRing"+i);
      }

      for(int i=1; i<=maxBondDist-1; i++)
      {  dNames.add("bdType"+i);
         dNames.add("bdRingSize"+i);
      }

      return dNames.toArray(new String[dNames.size()]);
   }
}

class PathDescription
{  private final OEMolBase mol;
   private final OEAtomBase[] pathAtoms;
   private final OEUnaryBondPred bondIsInRingFunctor = new OEBondIsInRing();

   String pathSmi = null;
   Map<String,Double> desc = new HashMap<String, Double>();

   public PathDescription(OEMolBase mol, OEAtomBase at1, OEAtomBase at2)
   {  pathAtoms = CATSIndexer.getAtoms(oechem.OEShortestPath(at1, at2));
      this.mol = mol;
   }

   public double getDescriptor(String name)
   {  if( desc.containsKey(name) )
         return desc.get(name);
      return 0D;
   }

   int getPathLength()
   {  return pathAtoms.length;
   }

   void computeDescription()
   {  double vdwSum = 0;
      double bondOrderSum = 0;
      int nBonds = 0;
      int nRingBonds = 0;
      int nSingle = 0;
      int nDouble = 0;
      int nCisBond = 0;
      int nTransBond = 0;
      int nTripple = 0;
      int nArom= 0;

      int totSubstitutions = 0;
      int totValence = 0;

      int atNum = 0;
      OEAtomBase lastAt = null;
      for(OEAtomBase at: pathAtoms)
      {  atNum++;

         vdwSum += Atom.VDW_RADIUS[at.GetAtomicNum()];

         int nSubstiution = at.GetExplicitDegree();
         int valence = at.GetValence();
         totSubstitutions += nSubstiution;
         totValence += valence;

         addDesc(desc,"atVal",    atNum,valence);
         addDesc(desc,"atSubst",  atNum,nSubstiution);
         addDesc(desc,"atIsRing", atNum, at.IsInRing() ? 1 : 0);

         if(lastAt != null)
         {  nBonds++;
            OEBondBase bd = at.GetBond(lastAt);
            int bdType = 4;
            if( ! bd.IsAromatic() ) bdType = bd.GetIntType();
            bondOrderSum += (bdType == 4 ? 1.5 : bdType);
            addDesc(desc, "bdType", nBonds, bdType == 4 ? 1.5 : bdType);

            if( bdType == 1 ) nSingle++;
            else if( bdType == 2 ) nDouble++;
            else if( bdType == 3 ) nTripple++;
            else if( bdType == 4 ) nArom++;

            if( bd.IsInRing() ) nRingBonds++;

            addDesc(desc,"bdRingSize", nBonds, oechem.OEBondGetSmallestRingSize(bd));

            if( atNum > 2 && pathAtoms.length >= atNum + 1)
            {  if( bd.HasStereoSpecified() )
               {  OEAtomBaseVector v = new OEAtomBaseVector();
                  v.add(pathAtoms[atNum-3]);
                  v.add(pathAtoms[atNum]);
                  boolean isTrans =
                        bd.GetStereo(v, OEBondStereo.CisTrans) == OEBondStereo.Trans;
                  if( isTrans )
                     nTransBond++;
                  else
                     nCisBond++;

               } else if( isTransAmide(bd, pathAtoms[atNum],pathAtoms[atNum-1], pathAtoms[atNum-2], pathAtoms[atNum-3]))
               {  nTransBond++;

               } else if( bd.IsInRing() )
               {  double cisTrans = getRingBondGeom(bd, pathAtoms[atNum],pathAtoms[atNum-1], pathAtoms[atNum-2], pathAtoms[atNum-3]);
                  if( cisTrans < 0 )
                     nCisBond += cisTrans;
                  else
                     nTransBond += cisTrans;
               }
            }

         }
         lastAt = at;
      }
      desc.put("pathLen",      (double)pathAtoms.length);
      desc.put("vdwSum",       vdwSum);
      desc.put("bondOrderSum", bondOrderSum);
      desc.put("totVal",   (double)totValence);
      desc.put("totSubst", (double)totSubstitutions);

      desc.put("nRingBonds", (double)nRingBonds);
      desc.put("nSingle",    (double)nSingle);
      desc.put("nDouble",    (double)nDouble);
      desc.put("nCisBond",   (double)nCisBond);
      desc.put("nTransBond", (double)nTransBond);
      desc.put("nTripple",   (double)nTripple);
      desc.put("nArom",      (double)nArom);
   }


   /**
    *
    * @param bd   central bond
    * @return -1 if cis 0 if staggered/unknown, 1 if trans
    */
   private double getRingBondGeom(OEBondBase bd, OEAtomBase at1, OEAtomBase at2,
                                                  OEAtomBase at3, OEAtomBase at4)
   {  if( ! bd.IsInRing() ) return 0;

      // to make this really clean we would need to be able to identify if atoms
      // are in the same small ring. OE does not provide this info so we are
      // using approxiamtions

      OEBondBase bd1 = at1.GetBond(at2);
      OEBondBase bd3 = at3.GetBond(at4);
      boolean bd1IsRing = bd1.IsInRing();
      boolean bd3IsRing = bd3.IsInRing();
      int at3NRingBond = getNRingBonds(at3);

      if( ! bd1IsRing ) // bd1 is exoCyclic
      {  return processExoBondNextToBridge(bd, bd3, bd3IsRing, at3, at4, at3NRingBond);
      }

      int at2NRingBond = getNRingBonds(at2);
      if( ! bd3IsRing ) // bd2 is exoCyclic
      {  return processExoBondNextToBridge(bd, bd1, bd1IsRing, at2, at1, at2NRingBond);
      }

      // both bd1 and bd3 are in rings
      if(at2NRingBond == 2 && at3NRingBond == 2) return -1; // simple ring -> cis

      if(at2NRingBond == 2 && at3NRingBond == 3)   // simple ring to bridge head: [Du]12CC[Ag][B]C1[Du]CC2
      {  if(getNRingBonds(at4) == 3 ) return -1; // simple ring to annulated bond
         OEAtomBase otherRingAtom = getOtherRingAtom(at3, bd, bd3);
         if( getNRingBonds(otherRingAtom) == 3 ) return 1; // simple ring to exo-ring bond on bridgehead: C12CC[Ag][B]C1[Du]CC2

         return 0; // this is a more complex ring that we do not currently handle
      }

      if(at3NRingBond == 2 && at2NRingBond == 3)   // simple ring to bridge head: [Du]12CC[Ag][B]C1[Du]CC2
      {  if(getNRingBonds(at1) == 3 ) return -1; // simple ring to annulated bond
         OEAtomBase otherRingAtom = getOtherRingAtom(at2, bd, bd1);
         if( getNRingBonds(otherRingAtom) == 3 ) return 1; // simple ring to exo-ring bond on bridgehead: C12CC[Ag][B]C1[Du]CC2

         return 0; // this is a more complex ring that we do not currently handle
      }

      if(at2NRingBond == 3 && at3NRingBond == 3)  // bridge on annulated ring
      {  return 0; // this is an anlatead ring bridge head, currently we can not
                   // identify cis trans on this with the oechem toolkit.
      }

      return 0; // more complex situation e do not handle
   }

   /**
    * @param bd central bond
    * @param bd3 exit bond
    * @param bd3IsRing
    * @return -1 if cis, 1 if trans and 0 if not known or staggered
    */
   private double processExoBondNextToBridge(OEBondBase bd, OEBondBase bd3,
            boolean bd3IsRing, OEAtomBase at3, OEAtomBase at4, int at3NRingBond)
   {
      if( ! bd3IsRing ) return -1; // two exocyclic -> cis: [Ag]B1CCC[C@@H]1[Du]
         if( at3NRingBond == 2 ) return 1; // bd1 exo, bd3 in simple ring -> trans: [Ag]B1CC[Du][C@@H]1C
         if( at3NRingBond == 4 ) return 0; // bd1 exo, at3 spiro: [Ag]B1CC[Du]C12[Du][Du]2

         // bd1 exo, at3 is bridge head atom
         assert at3NRingBond == 3 : "Very strange structure: atom with > 4 ring bonds";

         int at4NRingBond = getNRingBonds(at4);
         if(at4NRingBond == 3) return 1; // bd1 is exo, bd3 is annulated bond -> trans: [Ag]B1CC[Du]2[C@@H]1CCC2

         OEAtomBase otherRingAt = getOtherRingAtom(at3, bd, bd3);
         if(getNRingBonds(otherRingAt) == 3) return -1; //bd1 is exo, bd3 is also exo but also in ring: [Ag]B1CCC2[C@@H]1[Du]CCC2

         return 0; // These are complex bridged systems that we can not percieve
   }

   /**
    * @return The neighboring ring atom on at that is neither attached by bd1 nore by bd2.
    */
   private OEAtomBase getOtherRingAtom(OEAtomBase at, OEBondBase bd1, OEBondBase bd2)
   {  OEBondBase otherRingBond = null;
      OEBondBaseIter bdIt = at.GetBonds(bondIsInRingFunctor);
      while(bdIt.hasNext())
      {  otherRingBond = bdIt.next();
         if( otherRingBond.GetIdx() != bd1.GetIdx() && otherRingBond.GetIdx() != bd2.GetIdx() )
            break;
      }
      bdIt.delete();
      return otherRingBond.GetNbr(at);
   }

   private static int getNRingBonds(OEAtomBase at)
   {  OEBondBaseIter bIt = at.GetBonds();
      int nRing = 0;
      while(bIt.hasNext())
         if( bIt.next().IsInRing() ) nRing++;
      bIt.delete();
      return nRing;
   }

   private static boolean isTransAmide(OEBondBase bd, OEAtomBase at1, OEAtomBase at2,
                                                      OEAtomBase at3, OEAtomBase at4)
   {  if( bd.IsInRing() ) return false;
      if( bd.GetIntType() != 1 )return false;
      if( at2.GetBond(at1).GetIntType() != 1) return false;
      if( at3.GetBond(at4).GetIntType() != 1) return false;

      if(at2.GetAtomicNum() == OEElemNo.N && at2.GetImplicitHCount() == 1)
      {  if(at3.GetAtomicNum() == OEElemNo.C)
            return CATSIndexer.CARBONYLSS.AtomMatch(at3);

      } else if(at3.GetAtomicNum() == OEElemNo.N && at3.GetImplicitHCount() == 1)
      {  if(at2.GetAtomicNum() == OEElemNo.C)
            return CATSIndexer.CARBONYLSS.AtomMatch(at2);
      }
      return false;
   }

   String getSmiles()
   {  if( pathSmi == null )
      {  OEBitVector bv = new OEBitVector(mol.GetMaxAtomIdx());
         for(OEAtomBase at: pathAtoms)
            bv.SetBitOn(at.GetIdx());
         OEAtomIdxSelected pathFunctor = new OEAtomIdxSelected(bv);
         OEGraphMol pathMol = new OEGraphMol();
         oechem.OESubsetMol(pathMol, mol, pathFunctor);

         pathSmi = oechem.OECreateSmiString(pathMol);

         pathMol.delete();
         pathFunctor.delete();
         bv.delete();
      }
      return pathSmi;
   }

   /**
    * Add a descriptor to the desc map for an object (atom/bond) with teh given id.
    * @param val value of descriptor
    */
   private static void addDesc(Map<String, Double> desc, String name, int id,
            double val)
   {  desc.put(name+id, val);
   }

   public void close()
   {  bondIsInRingFunctor.delete();
   }
}
