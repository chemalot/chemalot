%mem=#mem#
%nprocshared=1
%chk=#FName#.chk
#T freq HF/3-21g

#FName#

#XYZ#

--Link1--
%mem=#mem#
%nprocshared=1
%chk=#FName#.chk
#T Geom=AllCheck Guess=TCheck opt(readfc,ModRedundant,MaxCycles=1000)  HF/6-31g* NoSym

#FIX#

--Link1--
%mem=#mem#
%nprocshared=1
%chk=#FName#.chk
# Geom=AllCheck Guess=TCheck MP2/6-311+g**

#FIX#

