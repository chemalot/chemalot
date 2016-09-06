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

package com.gNova.circularFP;

import java.util.*;
import java.util.Map.Entry;

import openeye.oechem.*;

/**
 *
 * @author tjodonnell
 *
 */
public class CFP
{  private final boolean debug;

   private final HashMap<Integer, Integer> fp;
   private final List<BondIdxPerAtom> bondIdxListsPerLevel = new ArrayList<BondIdxPerAtom>();

   int maxAtomIdx;

   private String[] stdSmarts = {
   "[a] aromatic",
   "[F,Cl,Br,I] halogen",
   "[$([$([#8,#16]);!$(*=NO);!$(*N=O);X1,X2]),$([#7;v3;!$([nX3]);!$(*(-a)-a)])] acceptor",
   "[$([O;H1,-&!$(*-N=O)]),$([S;H1&X2,-&X1]),$([#7;!H0;!$(*(S(=O)=O)C(F)(F)F);!$(n1nnnc1);!$(n1nncn1)]),$([#7;-])] donor",
   "[$([O;H1]-[C,S,P]=O),$([*;-;!$(*[*;+])]),$([N;!H0](S(=O)=O)C(F)(F)F),$(n1nnnc1),$(n1nncn1)] acidic",
   "[$([NH2]-[CX4]),$([NH](-[CX4])-[CX4]),$(N(-[CX4])(-[CX4])-[CX4]),$([*;+;!$(*[*;-])])$(N=C-N),$(N-C=N)] basic"
   };

   private List<String>Smarts;
   private OESubSearch[] functionalQuery;

   public CFP(boolean debug)
   {  this.debug = debug;
      this.fp  = new HashMap<Integer,Integer>();
   }

   public void initializeSmarts(List<String>smarts)
   {  if (smarts != null && smarts.size() > 0)
      {  this.Smarts = smarts;
      } else
      {  this.Smarts = new ArrayList<String>(Arrays.asList(stdSmarts));
      }
      functionalQuery = new OESubSearch[this.Smarts.size()];
        int i=0;
        for (String s : this.Smarts) {
           this.functionalQuery[i] = new OESubSearch(s);
           //this.functionalQuery[i].SetMaxMatches(1);
           ++i;
      }
      if (debug) for (String s: Smarts) System.err.println(s);
   }

   public void clear()
   {  this.fp.clear();
   }

/*
   private int hash(List<Integer>alist) {
   return alist.hashCode();
   }
*/
   private int hash(List<Integer> value)
   {
      if (value == null)
         return 0;

      long aprime = 2868947;
      long result = 0;
      for (int v : value)
      {  result = result*aprime + v;
      }
      return (int)result;
  }

   private List<Integer> initialFunctionalIdentifiers(OEGraphMol mol)
   {  int[] aid =  new int[maxAtomIdx+1];
      int imatch = 0;
      //this.initializeSmarts(null);
      for (OESubSearch ss : functionalQuery)
      {  if (this.debug)
         {  String s = this.Smarts.get(imatch);
            int i = s.indexOf(" ");
            if (i == -1)
            {  System.err.print(s + ":");

            } else
            {  System.err.print(s.substring(i+1) + ":");
            }
         }

         OEMatchBaseIter ssIt = ss.Match(mol);
         while( ssIt.hasNext() )
         {  OEMatchBase match = ssIt.next();

            OEMatchPairAtomIter aPairIt = match.GetAtoms();
            while( aPairIt.hasNext() )
            {  OEMatchPairAtom ma = aPairIt.next();

               int iatom = ma.getTarget().GetIdx();
               if (this.debug) System.err.print(iatom + " ");
               aid[iatom] += 1<<imatch;
            }
            aPairIt.delete();
         }
         ssIt.delete();

         if (this.debug) System.err.println(" ");
         ++imatch;
      }

      List<Integer> atomid = new ArrayList<Integer>(maxAtomIdx+1);
      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase atom = aIt.next();
         if (atom.IsHydrogen())
            atomid.add(null);
         else
            atomid.add(aid[atom.GetIdx()]);
      }
      aIt.delete();

      return atomid;
   }

   private List<Integer> initialAtomicIdentifiers(OEGraphMol mol)
   {
      List<Integer> atomid = new ArrayList<Integer>(maxAtomIdx+1);
      List<Integer> afp = new ArrayList<Integer>();
      int i = 0;
      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase atom = aIt.next();

         afp.clear();
         if (atom.IsHydrogen())
         {  atomid.add(null);

         } else
         {  int hcount = (atom.GetExplicitHCount() + atom.GetImplicitHCount());
            int atnum = atom.GetAtomicNum();
            afp.add((atom.IsInRing()   ? 1 : 0));
            afp.add((atom.IsAromatic() ? 1 : 0));
            //afp.add(atom.GetExplicitDegree());
            afp.add(atom.GetHvyDegree());
            afp.add(atom.GetHvyValence());
            afp.add(hcount);
            afp.add(atom.GetFormalCharge());
            afp.add(atom.GetAtomicNum());
            afp.add(Math.round(oechem.OEGetDefaultMass(atnum)));
   /*
            int atomID;
            atomID = (((((((
                  (atom.IsInRing() ? 1 : 0)
               *  2) + (atom.IsAromatic() ? 1 : 0)
               *  5) + atom.GetExplicitDegree()
               *  5) + atom.GetHvyValence()
               *  4) + hcount
               *  5) + (atom.GetFormalCharge()+2)
               * 54) + atom.GetAtomicNum()
               *128) + Math.round(oechem.OEGetDefaultMass(atnum));
            atomid.add(atomID);
            if (this.debug) System.err.println(atomID);
   */
            if (this.debug) System.err.println(++i+" "+afp+" "+this.hash(afp));
            atomid.add(this.hash(afp));
         }
      }
      aIt.delete();

      return atomid;
   }

   /**
    * Generate the FP0 initial identifier for each atom
    * @param mol
    * @param type
    * @return list of atomIDs for derived only from atom and its bonds, one id
    *         per atom index by atomIndex.
    */
   private List<Integer> initialAtomIdentifiers(OEGraphMol mol, String type)
   {  List<Integer> atomid;
      if (type.equals("functional"))
      {  atomid = initialFunctionalIdentifiers(mol);

      } else
      {  atomid = initialAtomicIdentifiers(mol);
      }
      return atomid;
   }

   /** Compute atom identifiers for this iteration (iter).
    * This follows description on page 744 of Rogers & Hahn paper.
    * @param mol molecule object on which to compute next iterations atom identifiers.
    * @param oldAtIds list of atom identifiers of previous iteration.
    * @param level cfp iteration number
    * @return list of atom identifier for this iteration
    */
   private List<Integer> iterativeAtomIdentifiers(OEGraphMol mol, List<Integer> oldAtIds, int level)
   {  List<Integer> newAtIds = new ArrayList<Integer>(maxAtomIdx+1);

      /** list of neighbor bondOrder AtIds to be hashed into new atID for current atom*/
      List<Integer> oid = new ArrayList<Integer>();

      /** bondid list of atomIdx by bond orders 0,1,2,4(aromatic) */
      int single = 0;
      int dble = 1;
      int tripple = 2;
      int arom = 3;

      List<List<Integer>> bondid = new ArrayList<List<Integer>>();
      bondid.add(new ArrayList<Integer>());
      bondid.add(new ArrayList<Integer>());
      bondid.add(new ArrayList<Integer>());
      bondid.add(new ArrayList<Integer>());
      int iatom = 0;

      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase atom = aIt.next();

         if (atom.IsHydrogen())
         {  newAtIds.add(null);
            continue;
         }

         oid.clear();
         oid.add(level);
         oid.add(oldAtIds.get(atom.GetIdx()));

         bondid.get(single).clear();
         bondid.get(dble).clear();
         bondid.get(tripple).clear();
         bondid.get(arom).clear();

         OEBondBaseIter bit = atom.GetBonds();
         while( bit.hasNext() )
         {  OEBondBase bd = bit.next();
            OEAtomBase nbr = bd.GetNbr(atom);
            if (!nbr.IsHydrogen())
            {  int bo = bd.IsAromatic() ? 4 : bd.GetOrder();
               bondid.get(bo-1).add(oldAtIds.get(nbr.GetIdx()));
            }
         }
         bit.delete();

         // compute list of neighbors atomIDs by bond orders
         for (int bdType=0; bdType<4; ++bdType)
         {  if (bondid.get(bdType).size() > 0)
            {  Collections.sort(bondid.get(bdType));

               for (Integer xid : bondid.get(bdType) )
               {  oid.add(bdType+1);
                  oid.add(xid);
               }
            }
         }
         if (this.debug) System.err.println((iatom++)+" "+oid+" "+this.hash(oid));
         // has list of atomids and bond orders into next atomID
         newAtIds.add(hash(oid));
      }
      aIt.delete();

      return newAtIds;
   }

   /**
    * Check for fragments which are duplicates in terms of containing the same
    * bonds.
    * This is described on page 745 of the Rogers & Hahn paper.
    *
    * This could be possibly much simpler based purley on bondifexes??
    */
   private void addNonDuplicates(int level, List<List<Integer>> atomid)
   {  int nnew = 0;
      List<Integer> aid = atomid.get(level);
      List<Set<Integer>> curBondIdxByAtomIdx = bondIdxListsPerLevel.get(level).bondIdxPerAtomIdx;

      if (this.debug)
      {  System.err.println("Iteration "+level);
         System.err.println(curBondIdxByAtomIdx);
      }

      for(int atIdx=0; atIdx < curBondIdxByAtomIdx.size(); atIdx++)
      {  Set<Integer> atomBondSet = curBondIdxByAtomIdx.get(atIdx);
         if( atomBondSet == null ) continue; //hydrogen or null

         // check if bondset was present in previous level
         boolean isDuplicate = false;
         for( int prevLevel = level-1; prevLevel >=0; prevLevel-- )
         {  BondIdxPerAtom prevBondIdxByAtomIdx = bondIdxListsPerLevel.get(prevLevel);
            if( prevBondIdxByAtomIdx.contains(atomBondSet))
            {  if (debug) System.err.println("dup from previous iteration " + prevLevel);
               isDuplicate = true;
               break;
            }
         }

         if( isDuplicate ) continue; // duplicate due to previous level

         // check for duplicates in this iteration
         Integer minAtId = aid.get(atIdx);
         for(int atIdx2=0; atIdx2 < curBondIdxByAtomIdx.size(); atIdx2++)
         {  if( atIdx2 == atIdx ) continue;

            if( atomBondSet.equals( curBondIdxByAtomIdx.get(atIdx2)))
            {  Integer otherAtId = aid.get(atIdx2);

               if( minAtId.intValue() > otherAtId.intValue() )
               {  if (debug)
                     System.err.println("Identical ignored:"+ minAtId +" > than " + otherAtId );

                  isDuplicate = true;
                  break;

               } else if( minAtId.intValue() == otherAtId.intValue() )
               {  if( atIdx > atIdx2 )
                  {  if (debug)
                        System.err.println("Identical ignored:"+ minAtId +" == " + otherAtId );

                     isDuplicate = true;
                     break;
                  }

                  assert atIdx < atIdx2;

                  if( debug )
                     System.err.println("new with dupl "+ otherAtId +" == " + minAtId);

               } else if( debug && minAtId.intValue() < otherAtId.intValue() )
               {  System.err.println("new dupl "+ minAtId + " < " + otherAtId );
               }
            }
         }

         if( isDuplicate ) continue; // duplicate due to previous level

         Integer orgCount = fp.get(minAtId);
         if( orgCount == null )
            fp.put(minAtId, Integer.valueOf(1));
         else
            fp.put(minAtId, ++orgCount);

         ++nnew;
      }

      if (this.debug)
      {  HashSet<Integer>uniq = new HashSet<Integer>(aid);
         System.err.println(nnew + " new" + " " + uniq.size() + " unique");
         System.err.println("fp:"+this.fp.size()+" "+this.fp);
      }
   }

   public void generate(OEGraphMol mol, int level, String fptype)
   {  maxAtomIdx                  = mol.GetMaxAtomIdx();
      List<List<Integer>> atomid  = new ArrayList<List<Integer>>(level);
      atomid.add(initialAtomIdentifiers(mol, fptype));

      fp.clear();
      bondIdxListsPerLevel.clear();

      incrementCounts(atomid.get(0));
      if (this.debug)
      {  System.err.println("Iteration 0");
         System.err.println(" fp:"+this.fp);
      }

      // keep list of bonds of fragments at this level
      BondIdxPerAtom curBondListByAtomIdx = BondIdxPerAtom.createLevel0(mol, debug);
      bondIdxListsPerLevel.add(curBondListByAtomIdx);

      for (int iter=1; iter<=level; ++iter)
      {
         atomid.add( iterativeAtomIdentifiers(mol, atomid.get(iter-1), iter));

         curBondListByAtomIdx = curBondListByAtomIdx.createNextLevel();
         bondIdxListsPerLevel.add(curBondListByAtomIdx);

         if (iter > 0) addNonDuplicates(iter, atomid);
      }
   }

   private void incrementCounts(List<Integer> atIds)
   {  for( Integer atId : atIds)
      {  if( atId == null ) continue;  // hydrogen or deleted atom

         Integer orgCount = fp.get(atId);
         if( orgCount == null )
            fp.put(atId, Integer.valueOf(1));
         else
            fp.put(atId, ++orgCount);
      }
   }

   public Set<Integer> get()
   {  //return new HashSet<Integer>(this.fp);
      return fp.keySet();
   }

   /**
    * Get map of (fragmentID, countOccurence) for last generated fingerprint.
    * @param foldToSize if 0 no folding, if > 0 folded into [0, foldToSize-1].
    */
   public Map<Integer, Integer> getCounts(int foldToSize)
   {  if( foldToSize == 0 )
         return Collections.unmodifiableMap(fp);

      Map<Integer,Integer> folded = new HashMap<Integer, Integer>(fp.size());
      for( Entry<Integer, Integer> frag : fp.entrySet())
      {  int atId = frag.getKey();
         atId = Math.abs(atId) % foldToSize;
         Integer orgVal = folded.get(atId);
         if( orgVal != null )
            folded.put(atId, orgVal + frag.getValue());
         else
            folded.put(atId, frag.getValue());
      }
      return folded;
   }
}


/**
 * Helper class containing for each atom in a molecule and a given cfp level
 * all bonds for each fragment indexed by atomIdx.
 *
 * @author albertgo
 *
 */
class BondIdxPerAtom
{  /** bondIdxPerAtomIdx.get(n) =
    * list of sets of bond indexes included in fragment of with atomIdx = n.
    *
    * It also contains a list of terminal atoms per atomIdx centered fragment
    * so that in {@see #createNextLevel()) we can add the next level of bonds.
    */
   final List<Set<Integer>> bondIdxPerAtomIdx;

   /** terminalAtomsByAtIdx.get(atomIdx) :list of terminal atoms centered around
    * atomIdx generated by the breath first search at this level.
    */
   final List<Set<OEAtomBase>> terminalAtomsByAtIdx;
   final int level;
   private final int maxAtomIdx;
   private final boolean debug;

   private BondIdxPerAtom(int level, int maxAtomIdx, boolean debug)
   {  this.level = level;
      this.maxAtomIdx = maxAtomIdx;
      this.debug = debug;
      bondIdxPerAtomIdx = new ArrayList<Set<Integer>>(maxAtomIdx+1);
      terminalAtomsByAtIdx = new ArrayList<Set<OEAtomBase>>(maxAtomIdx+1);
   }

   /** used to create a deep copy of the BondIdxPerAtom object
    */
   private BondIdxPerAtom(int level, BondIdxPerAtom src )
   {  this.level = src.level;
      this.maxAtomIdx = src.maxAtomIdx;
      this.debug = src.debug;
      bondIdxPerAtomIdx = new ArrayList<Set<Integer>>(maxAtomIdx+1);
      terminalAtomsByAtIdx = new ArrayList<Set<OEAtomBase>>(maxAtomIdx+1);

      bondIdxPerAtomIdx.addAll( src.bondIdxPerAtomIdx );
      terminalAtomsByAtIdx.addAll( src.terminalAtomsByAtIdx);

      for( int atIdx = 0; atIdx < src.bondIdxPerAtomIdx.size(); atIdx++ )
      {  Set<Integer> srcBSet = src.bondIdxPerAtomIdx.get(atIdx);
         if(srcBSet == null ) continue;


         HashSet<Integer> newBSet = new HashSet<Integer>();
         newBSet.addAll(srcBSet);
         bondIdxPerAtomIdx.set(atIdx, newBSet);

         Set<OEAtomBase> srcTSet = src.terminalAtomsByAtIdx.get(atIdx);
         HashSet<OEAtomBase> newTSet = new HashSet<OEAtomBase>();
         newTSet.addAll(srcTSet);
         terminalAtomsByAtIdx.set(atIdx, newTSet);
      }
   }

   /**
    * Create empty list of bonds with atoms themselfs as terminal atoms.
    *
    */
   static BondIdxPerAtom createLevel0(OEGraphMol mol, boolean debug)
   {  BondIdxPerAtom newBondIdxPerAtomIdxList
         = new BondIdxPerAtom(0, mol.GetMaxAtomIdx(), debug);

      OEAtomBaseIter aIt = mol.GetAtoms();
      while( aIt.hasNext() )
      {  OEAtomBase at = aIt.next();
         if( at.IsHydrogen() ) continue;

         int atIdx = at.GetIdx();

         // ensure capacity
         while( newBondIdxPerAtomIdxList.bondIdxPerAtomIdx.size() <= atIdx )
         {  newBondIdxPerAtomIdxList.bondIdxPerAtomIdx.add(null);  // hydrogen or empty atomIdx
            newBondIdxPerAtomIdxList.terminalAtomsByAtIdx.add(null);
         }

         HashSet<Integer>    bSet   = new HashSet<Integer>(0); // no bonds for 0 level
         newBondIdxPerAtomIdxList.bondIdxPerAtomIdx.add(atIdx, bSet);

         HashSet<OEAtomBase> tAtoms = new HashSet<OEAtomBase>(1); // just itself
         tAtoms.add(at);
         newBondIdxPerAtomIdxList.terminalAtomsByAtIdx.add(atIdx,tAtoms);
      }
      aIt.delete();

      return newBondIdxPerAtomIdxList;
   }


   BondIdxPerAtom createNextLevel()
   {  BondIdxPerAtom newBondIdxPerAtom = new BondIdxPerAtom(level+1, this);
      List<Set<OEAtomBase>> newTerminalAtomsByAtIdx = newBondIdxPerAtom.terminalAtomsByAtIdx;
      List<Set<Integer>> newBondIdxPerAtomIdx = newBondIdxPerAtom.bondIdxPerAtomIdx;

      for( int atIdx=0; atIdx < terminalAtomsByAtIdx.size(); atIdx++)
      {  Set<OEAtomBase> tAtomSet = terminalAtomsByAtIdx.get(atIdx);
         if( tAtomSet == null ) continue; // null or hydrogen
         Set<OEAtomBase> newTAtomSet = new HashSet<OEAtomBase>();

         Set<Integer> bSet = newBondIdxPerAtomIdx.get(atIdx);

         for( OEAtomBase tAtom : tAtomSet )
         {  OEBondBaseIter bIt = tAtom.GetBonds();
            while( bIt.hasNext() )
            {  OEBondBase bd = bIt.next();
               Integer bIdx = bd.GetIdx();
               if( bSet.contains(bIdx) )  continue;

               // add new bond to bondlist for this atom
               bSet.add(bIdx);
               newTAtomSet.add(bd.GetNbr(tAtom));
            }
            bIt.delete();
         }

         // update terminalAtomList
         newTerminalAtomsByAtIdx.set(atIdx, newTAtomSet);
         if( debug )
         {  System.err.println("AtIdx "+ atIdx);
            System.err.println("\tbd  "+ newBondIdxPerAtom.bondIdxPerAtomIdx.get(atIdx));
            System.err.println("\ttat "+ newBondIdxPerAtom.terminalAtomsByAtIdx.get(atIdx));
         }
      }

      return newBondIdxPerAtom;
   }

   /*
    * True if bSet is present on any of the atoms in this BondIdxPerAtom
    */
   boolean contains(Set<Integer> bSet)
   {  for(Set<Integer> thisBSet : bondIdxPerAtomIdx )
         if(bSet.equals(thisBSet)) return true;

      return false;
   }
}
