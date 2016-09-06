# updateing an existing model in development:
1) edit "properties.xml" and change date enetity for model
2) edit model specific xml file and veryfy that descriptors are still the same
   change if needed
3) test
    setDev;echo CCO|sdfCalcProps.csh -in .smi -out .sdf cCYP_Inh
4) If everithing works proceed to "releasing changes to prd"
5) svn commit -m <"Commit message explaining changes maade">


# releasing changes to prd:

1) sudo -u smdi tcsh -l     #login as user smdi
2) make sure ~/dev/Aestel/config/properties contains the valid configuration files
  eg. vi properties.xml -> change date in ENTITY for ppb
2.1) test: setDev; echo CCO | sdfCalcProps.csh -in .smi -out .sdf cPPB
2.2) commit changes to dev svn: 
   cd ~smdi/dev/Aestel/config/properties
   # list all changes
   svn status
   # add any required new files whih is not in svn yet (marked with ?):
   svn add <filenames>
   # commit files with relevant changes:
   svn commit -m "<change message>" <filenames>
   # delete any junk files in the directory
3) cd ~/prd/Aestel
4) Maximize your putty window. To compare and copy development code from dev to prd execute:
 repCP.pl config/properties
  Hit C: to mark file for copy
  Hit I: to ignore file
  Hit V: to see file with vimdiff
  
     in diffs development is on left side and prd on right side
     vimdiff commands:
     ctrl-w ctrl-w    change window side
     do  use block from other side of window
     dp  put block from this window into other window
5) repCP.pl will have created file named repCPComms<date>.csh in ~/prd/Aestel
  vi file and isolate lines which you want to execute
  save those lines as t
  source t
6)  setPrd ; echo CCO | sdfCalcProps.csh -in .smi -out .sdf <ModelNames>
7) svn commit -m "prd release of first sdfRModel cLM and cHep models, congratulations everybody!" config/properties
