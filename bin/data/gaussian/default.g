%mem=#mem#
%nprocshared=4
%chk=#FName#.chk
# opt(MaxCycles=1000,tight) freq AM1 NoSym

#FName#

#XYZ#

--Link1--
%mem=#mem#
%nprocshared=2
%chk=#FName#.chk
# Geom=AllCheck Guess=TCheck opt(MaxCycles=500,readfc) HF/3-21g*  NoSym

--Link1--
%mem=#mem#
%nprocshared=2
%chk=#FName#.chk
# Geom=AllCheck Guess=TCheck opt(MaxCycles=500,readfc) HF/6-31+g** NoSym

--Link1--
%mem=#mem#
%nprocshared=2
%chk=#FName#.chk
%Mem=30GB
# Geom=AllCheck Guess=TCheck mp2/6-31+g** NoSym

--Link1--
%mem=#mem#
%nprocshared=2
%chk=#FName#.chk
# Geom=AllCheck Guess=TCheck opt(MaxCycles=500,readfc) HF/6-31+g** scrf=iefpcm NoSym

--Link1--
%mem=#mem#
%nprocshared=2
%chk=#FName#.chk
# Geom=AllCheck Guess=TCheck MP2/6-31+g** scrf=iefpcm NoSym

