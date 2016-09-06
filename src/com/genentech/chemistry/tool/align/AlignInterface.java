package com.genentech.chemistry.tool.align;

import openeye.oechem.OEMolBase;

public interface AlignInterface
{  public void align(OEMolBase fitmol);
   public void close();
}
