%mem=10GB
%nprocshared=3
%chk=id1.0216_001.chk
#T freq HF/3-21g

id1.0216_001

0 1
 C     1.1670 0.6924 -0.0383
 C     2.4135 1.4623 -0.0820
 O     2.7292 2.4706 0.4868
 O     1.1933 -0.5056 -0.1027
 H     1.3639 -0.8697 0.7691
 H     3.5812 2.7695 0.1605
 H     2.5358 1.7058 -1.1374
 H     3.0190 0.7760 0.5099
 H     0.5550 1.0382 -0.8714
 H     0.7971 0.8777 0.9701


--Link1--
%mem=10GB
%nprocshared=3
%chk=id1.0216_001.chk
#T Geom=AllCheck Guess=TCheck opt(readfc,ModRedundant,MaxCycles=1000)  HF/6-31g* NoSym

D 4 1 2 3  F


--Link1--
%mem=10GB
%nprocshared=3
%chk=id1.0216_001.chk
# Geom=AllCheck Guess=TCheck MP2/6-311+g**

D 4 1 2 3  F


