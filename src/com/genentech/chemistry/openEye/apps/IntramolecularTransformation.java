package com.genentech.chemistry.openEye.apps;

import java.io.File;
import java.io.IOException;
import java.util.*;

import com.genentech.oechem.tools.OETools;

import openeye.oechem.*;

/*
 * Define a class to store the reaction information
 */
public class IntramolecularTransformation {
   private String smirks;
   private String rxnFileName;
   private OEMolBase reactant = null;
   private OEMolBase product = null;
   private OEQMol reactantQMol = null;
   private OEQMol productQMol = null;

   public String getSmirks() {
      return this.smirks;
   }

   public String getRxnFile() {
      return this.rxnFileName;
   }

   public OEMolBase getReactant() {
      return this.reactant;
   }

   public OEMolBase getProduct() {
      return this.product;
   }

   public OEQMol getReactantQMol() {
      return this.reactantQMol;
   }

   public OEQMol getProductQMol() {
      return this.productQMol;
   }

   public void delete() {
      if (this.reactant != null)
         this.reactant.delete();

      if (this.product != null)
         this.product.delete();

      if (this.reactantQMol != null)
         this.reactantQMol.delete();

      if (this.productQMol != null)
         this.productQMol.delete();
   }

   public List<OEMolBase> run(List<OEMolBase> targets, boolean makeHExplicit, boolean matchStereo)
         throws Exception {
      OEQMol qmol = new OEQMol();

      oemolistream rfile = new oemolistream(this.getRxnFile());

      int opt = OEMDLQueryOpts.ReactionQuery;
      if( ! makeHExplicit )
         opt |= OEMDLQueryOpts.SuppressExplicitH;

      if (matchStereo)
         opt |= OEMDLQueryOpts.MatchAtomStereo;

      if (!oechem.OEReadMDLReactionQueryFile(rfile, qmol, opt)) {
         throw new Exception("Unable to read reaction file: "
               + this.getRxnFile());
      }

      OELibraryGen lg = new OELibraryGen();

      if (!lg.Init(qmol)) {
         throw new Exception("Failed to initialize reaction with "
               + this.getRxnFile());
      }

      if( makeHExplicit )
      {  lg.SetExplicitHydrogens(true);
         lg.SetValenceCorrection(false);
      } else
      {  lg.SetExplicitHydrogens(false);
         lg.SetValenceCorrection(true);
      }

      List<OEMolBase> products = new ArrayList<OEMolBase>();

      for (OEMolBase mol : targets) {
         lg.SetStartingMaterial(mol, 0, true);

         Set<String> prdSet = new HashSet<String>();
         for (OEMolBase product : lg.GetProducts()) {

            // only store transformation which yielded unique products for one input
            OEGraphMol prd = new OEGraphMol(correctValence(product));
            String smi = OETools.molToCanSmi(prd, true);

            if( prdSet.add(smi) )
               products.add(prd);
         }
      }
      lg.delete();
      qmol.delete();

/*
      if (lg != null)
         lg.delete();
*/
      /*
      OEUniMolecularRxn umr = new OEUniMolecularRxn();
      if (!umr.Init(qmol)) {
         throw new Exception("Failed to initialize reaction with "
               + this.getRxnFile());
      }

      List<OEMolBase> products = new ArrayList<OEMolBase>();

      for (OEMolBase mol : targets) {
         System.err.println(OETools.molToCanSmi(mol, true));
         if (umr.constCall(mol)) {
            System.err.println(OETools.molToCanSmi(mol, true));
            System.err.println(OETools.molToCanSmi(correctValence(mol), true));
            products.add(this.correctValence(mol));
         }
      }

      if (umr != null)
         umr.delete();

*/
      rfile.close();
      rfile.delete();

      return products;
   }

   public IntramolecularTransformation(OEMolBase reactant,
         OEMolBase product, String smirks, String rxnFileName) {
      this.reactant = reactant;
      this.product = product;
      this.smirks = smirks;
      this.rxnFileName = rxnFileName;
   }

   /**
    * Simple utility to read molecule
    *
    * @param fileName
    * @return
    */
   private static OEMolBase readMolecule(String fileName) {
      oemolistream ifs = new oemolistream(fileName);
      OEMolBase mol = new OEGraphMol();
      oechem.OEReadMolecule(ifs, mol);
      ifs.close();
      ifs.delete();

      return mol;
   }

   /**
    * Simple utility to write molecule
    *
    * @param fileName
    * @param mol
    */
   private static void writeMolecule(String fileName, OEMolBase mol) {
      oemolostream ofs = new oemolostream(fileName);
      oechem.OEWriteMolecule(ofs, mol);
      ofs.close();
      ofs.delete();
   }

   public IntramolecularTransformation(String rxnFileName, boolean makeHExplicit) throws IOException {
      this.rxnFileName = rxnFileName;

      oemolistream ifs = new oemolistream(rxnFileName);

      OEGraphMol mol = new OEGraphMol();
      oechem.OEReadRxnFile(ifs, mol);
      this.smirks = oechem.OECreateCanSmiString(mol);

      this.reactant = new OEGraphMol();
      oechem.OESubsetMol(reactant, mol, new OEAtomIsInReactant());
      OEAtomBaseIter ait = reactant.GetAtoms();
      while(ait.hasNext())
         ait.next().SetRxnRole(OERxnRole.None);
      ait.delete();
      reactant.SetRxn(false);

      this.product =  new OEGraphMol();
      oechem.OESubsetMol(product, mol, new OEAtomIsInProduct());
      ait = product.GetAtoms();
      while(ait.hasNext())
         ait.next().SetRxnRole(OERxnRole.None);
      ait.delete();
      product.SetRxn(false);

      // Read in the MDL Query File
      int opts = OEMDLQueryOpts.ReactionQuery;
      if( ! makeHExplicit ) opts |= OEMDLQueryOpts.SuppressExplicitH;

      ifs.rewind();

      mol.Clear();;
      oechem.OEReadMDLReactionQueryFile(ifs, mol);
      OEQMolBase qmol = new OEQMol();
      oechem.OEBuildMDLQueryExpressions(qmol, mol, opts);

      this.reactantQMol = new OEQMol();
      oechem.OESubsetMol(reactantQMol, qmol, new OEAtomIsInReactant());

      this.productQMol = new OEQMol();
      oechem.OESubsetMol(productQMol, qmol, new OEAtomIsInProduct());

      mol.delete();
      ifs.close();
      ifs.delete();
   }

   /**
    * Old method just here until new is fully tested
    * @param rxnFileName
    */
   public IntramolecularTransformation(String rxnFileName, String thisMethodIsOld) throws IOException {
      this.rxnFileName = rxnFileName;

      oemolistream ifs = new oemolistream(rxnFileName);

      OEGraphMol mol = new OEGraphMol();
      oechem.OEReadRxnFile(ifs, mol);
      this.smirks = oechem.OECreateCanSmiString(mol);

      ifs.rewind();

      // Read the reactant
      OEGraphMol rmol = new OEGraphMol();
      oechem.OEReadMolecule(ifs, rmol);
      File tempReactantFile = File.createTempFile("reaction", ".mol");
      writeMolecule(tempReactantFile.toString(), rmol);
      this.reactant = readMolecule(tempReactantFile.toString());

      // Read the product
      OEGraphMol pmol = new OEGraphMol();
      oechem.OEReadMolecule(ifs, pmol);
      File tempProductFile = File.createTempFile("product", ".mol");
      writeMolecule(tempProductFile.toString(), pmol);
      this.product = readMolecule(tempProductFile.toString());

      ifs.close();
      ifs.delete();

      // Read in the MDL Query File
      int opts = OEMDLQueryOpts.Default | OEMDLQueryOpts.SuppressExplicitH;
      this.reactantQMol = new OEQMol();
      oemolistream rqfile = new oemolistream(tempReactantFile.toString());
      oechem.OEReadMDLQueryFile(rqfile, this.reactantQMol, opts);

      rqfile.close();
      rqfile.delete();

      this.productQMol = new OEQMol();
      oemolistream pqfile = new oemolistream(tempProductFile.toString());
      oechem.OEReadMDLQueryFile(pqfile, this.productQMol, opts);

      pqfile.close();
      pqfile.delete();

      tempReactantFile.delete();
      tempProductFile.delete();
   }

   private OEMolBase correctValence(OEMolBase mol)
   {
      OEAtomBaseIter ai = mol.GetAtoms();

      while (ai.hasNext()) {
         OEAtomBase at = ai.next();
         if (!(at.IsNitrogen() || at.IsOxygen() || at.IsSulfur() || at.IsCarbon()))
            continue;

         if (at.IsNitrogen()) {
            correctValence(at, 3);
         }
         if (at.IsOxygen()) {
            correctValence(at, 2);
         }
         if (at.IsCarbon()) {
            correctValence(at, 4);
         }
         if (at.IsSulfur()) {
            int d = at.GetDegree();
            int v = at.GetValence();

            switch (d - v) {
            case 0: correctValence(at, 2); break;
            case 1: correctValence(at, 4); break;
            case 2: correctValence(at, 6); break;
            default: break;
            }

         }
      }

      ai.delete();

      return mol;
   }

   public void correctValence(OEAtomBase atom, int targetValence)
   {
      int formalCharge = atom.GetFormalCharge();
      int nvalence = atom.GetValence();
      int nexplicitH = atom.GetExplicitHCount();
      int nimplicitH =atom.GetImplicitHCount();
      int valence = nvalence - formalCharge;

      if (valence == targetValence)
         return;
      if (valence < targetValence) {
         atom.SetImplicitHCount(Math.abs(nimplicitH + targetValence - valence));
      }
      else {
         int correctedImplicitHcount = Math.max(0,  targetValence - valence + nimplicitH);
         atom.SetImplicitHCount(correctedImplicitHcount);
         valence = atom.GetExplicitValence() + correctedImplicitHcount;
         if (valence > targetValence && nexplicitH >= valence - targetValence) {
            RemoveExplicitHydrogen(atom, valence - targetValence);
         }
      }
   }


   private static void RemoveExplicitHydrogen(OEAtomBase atom, int nremove)
   {
      OEMolBase mol = atom.GetParent();
      OEAtomBaseIter aiter = atom.GetAtoms();
      while (aiter.hasNext()) {
         OEAtomBase at = aiter.next();
         if (at.IsHydrogen()) {
            mol.DeleteAtom(at);
            nremove--;
            if (nremove == 0)
               break;
         }
      }
      aiter.delete();
   }

   public static void main(String args[]) {
      if (args.length < 3) {
         System.err.println("usage: IntramolecularTransformation <rxnFile> <input.mol> <output.mol>");
         System.exit(1);
      }

      try {
         IntramolecularTransformation it = new IntramolecularTransformation(
               args[0], false);

         List<OEMolBase> mols = new ArrayList<OEMolBase>();
         OEMolBase mol = readMolecule(args[1]);
         mols.add(mol);

         List<OEMolBase> tmols = it.run(mols, false, true);

         if (tmols.size() > 0) {
            writeMolecule(args[2], tmols.get(0));
         } else {
            System.err.println("No transformation");
         }
      } catch (Exception e) {
         System.err.println(e.getMessage());
         System.exit(1);
      }
   }
}
