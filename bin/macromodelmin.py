__doc__ = """

macromdelmin.py is a convenience script that wraps Schrodinger's macromodel bmin
executable.

For each structure in the input file, generate a MacroModel com file and run
the job.

Copyright Schrodinger, LLC. All rights reserved.

macromodelmin.py was written by Chris Higgs of Schrodinger, LLC.
and was added to chemalot by permission of Chris Higgs and Schrodinger, LLC.
on 5/13/2016

"""

import argparse
import os
import sys
import itertools
import math
import glob

#BEN
import subprocess

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.utils import fileutils, cmdline
from schrodinger.job import jobcontrol, queue


###############################################################################
# Globals and constants

FF_DICT = {'MM2*': 1, 'MM3*': 2, 'AMBER': 3, 'AMBER94': 4, 'OPLS': 5,
           'MMFF': 10, 'MMFFs': 10, 'OPLS_2001': 11, 'OPLS_2005': 14, 'OPLS_21':16}

MM_DICT = {'PRCG': 1, 'TNCG': 9, 'OSVM': 3, 'SD': 0, 'FMNR': 4, 'LBFGS': 10}

# BEN
BMIN = os.path.join(os.environ['SCHRODINGER'], 'bmin')


###############################################################################
def parse_args():
    """
    Parse the command line options.

    @return:  All script arguments and options.
    @rtype:  class:`argparse.Namespace`
    """

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=
                                     argparse.RawDescriptionHelpFormatter)

    parser.add_argument("infile",
                        help="Input file.")

    parser.add_argument("outfile",
                        help="Output file.")

    parser.add_argument("-forcefield",
                        type=str,
                        choices=['MM2*', 'MM3*', 'AMBER*', 'AMBER94', 'OPLS',
                                 'MMFF', 'MMFFs', 'OPLS_2001', 'OPLS_2005', 'OPLS_21'],
                        default="OPLS_2005",
                        help="ForceField. (Default=OPLS_2005).")

    parser.add_argument("-solvent",
                        type=str,
                        choices=['None', 'Water', 'Octanol', 'CHCL3'],
                        default="Water",
                        help="Solvent. (Default= Water).")

    parser.add_argument("-mini_method",
                        type=str,
                        choices=['PRCG', 'TNCG', 'OSVM', 'SD', 'FMNR', 'LBFGS'],
                        default="PRCG",
                        help="Minimization method. (Default=PRCG).")

    parser.add_argument("-mini_steps",
                        default="5000",
                        help="Minimization steps. (Default=5000).")

    parser.add_argument("-convergence",
                        default="0.0500",
                        help="Gradient convergence threshold. ("
                             "Default=0.0500).")

    parser.add_argument("-frozenHeavy",
                        action="store_true",
                        default=False,
                        help="Optimize only the hydrogen atoms. All heavy "
                             "atoms are frozen.")

    parser.add_argument("-frozenSmarts",
                        default="",
                        help="All atoms matching the SMARTS pattern are "
                             "frozen.")

    parser.add_argument("-fixedHeavy",
                        action="store_true",
                        default=False,
                        help="Optimize only the hydrogen atoms. All heavy "
                             "atoms are fixed using the specified force "
                             "constant (-constraintForce) and flat-bottom "
                             "width (-constraintWidth).")

    parser.add_argument("-fixedSmarts",
                        default="",
                        help="All atoms matching the SMARTS pattern are fixed "
                             "using the specified force constant "
                             "(-constraintForce) and flat-bottom width "
                             "(-constraintWidth).")

    parser.add_argument("-constrain_by_coord",
                        default="",
                        help="List of x,y,z coordinates for the four torsion "
                             "atoms. Must be in the following format "
                             "\"x,y,z; x,y,z; x,y,z; x,y,z\"")

    parser.add_argument("-constrain_by_index",
                        default="",
                        help="List of indices for the four torsion "
                             "atoms. Must be in a comma-separated list.")

    parser.add_argument("-constraintForce",
                        default="100.0",
                        help="Force constant for the constraint. ("
                             "Default=100.0).")

    parser.add_argument("-constraintWidth",
                        default="0.0",
                        help="Width of the flat-bottomed region, "
                             "in angstroms. Default=0.0).")

    parser.add_argument("-host",
                        metavar='<hostlist>',
                        default="localhost:1",
                        help="Run job remotely on indicated host entries. "
                             "<hostlist> is a comma-separated list of one or "
                             "more schrodinger host file entry names. The "
                             "number of processors the host may use for the "
                             "job can be specified after a colon. The default"
                             " number of processors is one. Example: "
                             "'name1:nprocs1,name2:nprocs2' "
                             "Default:localhost:1")

    args = parser.parse_args()

    return args


################################################################################
def check_args(cmd_args):
    """
    Check command-line options.

    @param cmd_args:  All script arguments and options.
    @type cmd_args:  class:`argparse.Namespace`

    @return:  True if constraints are used.
    @rtype:  bool
    """

    enable_constraints = False

    if not os.path.isfile(cmd_args.infile):
        print "Input file not found. Please ensure the file exists."
        sys.exit(1)

    if fileutils.get_structure_file_format(cmd_args.infile) != fileutils.SD:
        print "Error: Input file is not in SD format."
        sys.exit(1)

    if os.path.isfile(cmd_args.outfile):
        print "Output file exists. Please specify a different file name."
        sys.exit(1)

    if fileutils.get_structure_file_format(cmd_args.infile) != fileutils.SD:
        print "Error: Output file doesn't appear to be SD format."
        sys.exit(1)

    if cmd_args.frozenHeavy is True and len(cmd_args.frozenSmarts):
        print "Error: Cannot use both the -frozenHeavy and -frozenSmarts " \
              "options."
        sys.exit(1)

    if cmd_args.frozenHeavy is True and cmd_args.fixedHeavy is True:
        print "Error: Cannot use both the -frozenHeavy and -fixedHeavy options."
        sys.exit(1)

    if cmd_args.frozenHeavy is True and len(cmd_args.fixedSmarts):
        print "Error: Cannot use both the -frozenHeavy and -fixedSmarts " \
              "options."
        sys.exit(1)

    if cmd_args.frozenHeavy is True and len(cmd_args.constrain_by_coord):
        print "Error: Cannot use both the -frozenHeavy options and " \
              "-constrain_by_coord options."
        sys.exit(1)

    if cmd_args.frozenHeavy is True and len(cmd_args.constrain_by_index):
        print "Error: Cannot use both the -frozenHeavy options and " \
              "-constrain_by_coord options."
        sys.exit(1)

    if len(cmd_args.frozenSmarts) and cmd_args.fixedHeavy is True:
        print "Error: Cannot use both the -frozenSmarts and -fixedHeavy " \
              "options."
        sys.exit(1)

    if len(cmd_args.frozenSmarts) and len(cmd_args.fixedSmarts):
        print "Error: Cannot use both the -frozenSmarts and -fixedSmarts " \
              "options."
        sys.exit(1)

    if len(cmd_args.frozenSmarts) and len(cmd_args.constrain_by_coord):
        print "Error: Cannot use both the -frozenSmarts options and " \
              "-constrain_by_coord options."
        sys.exit(1)

    if len(cmd_args.frozenSmarts) and len(cmd_args.constrain_by_index):
        print "Error: Cannot use both the -frozenSmarts options and " \
              "-constrain_by_index options."
        sys.exit(1)

    if cmd_args.fixedHeavy is True and len(cmd_args.fixedSmarts):
        print "Error: Cannot use both the -fixedHeavy and -fixedSmarts and " \
              "options."
        sys.exit(1)

    if cmd_args.fixedHeavy is True and len(cmd_args.constrain_by_coord):
        print "Error: Cannot use both the -fixedHeavy and -constrain_by_coord" \
              " options."
        sys.exit(1)

    if cmd_args.fixedHeavy is True and len(cmd_args.constrain_by_index):
        print "Error: Cannot use both the -fixedHeavy and -constrain_by_index" \
              "options."
        sys.exit(1)

    if len(cmd_args.fixedSmarts) and len(cmd_args.constrain_by_coord):
        print "Error: Cannot use both the -fixedSmarts and " \
              "-constrain_by_coord options."
        sys.exit(1)

    if len(cmd_args.fixedSmarts) and len(cmd_args.constrain_by_index):
        print "Error: Cannot use both the -fixedSmarts and " \
              "-constrain_by_index options."
        sys.exit(1)

    if len(cmd_args.constrain_by_coord) and len(cmd_args.constrain_by_index):
        print "Error: Cannot use both the -constrain_by_coord and " \
              "-constrain_by_index options."
        sys.exit(1)

    if (cmd_args.frozenHeavy is True or len(cmd_args.frozenSmarts) or
            cmd_args.fixedHeavy is True or len(cmd_args.fixedSmarts) or
            len(cmd_args.constrain_by_coord) or
            len(cmd_args.constrain_by_index)):
            enable_constraints = True

    return enable_constraints


################################################################################
def create_com_file(inFilename, cmd_args, constraint_at_list):
    """
    Create the MacroModel com file.

    @param inFilename:  filename.
    @type inFilename:  string

    @param cmd_args:  All script arguments and options.
    @type cmd_args:  class:`argparse.Namespace`

   
    @param constraint_at_list:  List of atom indices to constrain of fix.
    @type constraint_at_list:  list

    @return:  Filenames
    @rtype:  str, str, str, str, str
    """

    line3 = None
    line4 = None
    line5 = None

    basename = fileutils.get_basename(cmd_args.outfile)

    comfile = basename + ".com"
    inmaefile = basename + "-in.mae"
    outmaefile = basename + "-out.mae"
    logfile = basename + ".log"
    
    writer = structure.StructureWriter(inmaefile)   
    for i, st in enumerate(structure.StructureReader(inFilename)):
        writer.append(st)
    writer.close()

    line1 = (" MMOD       0      1      0      0     0.0000     0.0000     "
             "0.0000     0.0000" + "\n")

    if cmd_args.forcefield != "MMFFs":
        line2 = (" FFLD" +
                 str(FF_DICT[cmd_args.forcefield]).rjust(8) +
                 "      1      0      0     1.0000     0.0000     0.0000     "
                 "0.0000" + "\n")
    else:
        line2 = (" FFLD" +
                 str(FF_DICT[cmd_args.forcefield]).rjust(8) +
                 "      1      0      1     1.0000     0.0000     0.0000     "
                 "0.0000" + "\n")

    if cmd_args.solvent == "None":
        line5 = (" BDCO       0      0      0      0    41.5692 99999.0000     "
                 "0.0000     0.0000" + "\n")

    elif cmd_args.solvent == "Water":
        line3 = (" SOLV       3      1      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
        line4 = (" EXNB       0      0      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
        line5 = (" BDCO       0      0      0      0    89.4427 99999.0000     "
                 "0.0000     0.0000" + "\n")

    elif cmd_args.solvent == "Octanol":
        line3 = (" SOLV       3      9      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
        line4 = (" EXNB       0      0      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
        line5 = (" BDCO       0      0      0      0    89.4427 99999.0000     "
                 "0.0000     0.0000" + "\n")

    elif cmd_args.solvent == "CHCL3":
        line3 = (" SOLV       3      5      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
        line4 = (" EXNB       0      0      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
        line5 = (" BDCO       0      0      0      0    89.4427 99999.0000     "
                 "0.0000     0.0000" + "\n")
    line5b    = (" BGIN       0      0      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
                 
    line6 = (" READ       0      0      0      0     0.0000     0.0000     "
             "0.0000     0.0000" + "\n")

    line7 = (" CONV       2      0      0      0" +
             str("%.4f" % float(cmd_args.convergence)).rjust(11) +
             "     0.0000     0.0000     0.0000" + "\n")

    if cmd_args.mini_method != "TNCG":
        line8 = (" MINI" + str(MM_DICT[cmd_args.mini_method]).rjust(8) +
                 "      0" + str(cmd_args.mini_steps).rjust(7) +
                 "      0     0.0000     0.0000     0.0000     0.0000" + "\n")
    else:
        line8 = (" MINI       9      1" + str(cmd_args.mini_steps).rjust(7) +
                 "      0     0.0000     0.0000     0.0000     0.0000" + "\n")

    line9     = (" END        0      0      0      0     0.0000     0.0000     "
                 "0.0000     0.0000" + "\n")
                 
    with open(comfile, "w") as com_file:
        com_file.write(inmaefile + "\n")
        com_file.write(outmaefile + "\n")
        com_file.write(line1)
        com_file.write(line2)

        if line3 is not None:
            com_file.write(line3)
            com_file.write(line4)
            com_file.write(line5)
        else:
            com_file.write(line5)

        com_file.write(line5b)
        com_file.write(line6)

        if len(constraint_at_list):
            for at in constraint_at_list:
                if (cmd_args.frozenHeavy is True or
                        len(cmd_args.frozenSmarts) or
                        len(cmd_args.constrain_by_coord)):
                    line = (" FXAT "'{:>7}'.format(at) +
                            "      0      0      0    -1.0000     0.0000     "
                            "0.0000     0.0000" + "\n")
                    com_file.write(line)

                elif cmd_args.fixedHeavy is True or len(cmd_args.fixedSmarts):
                    line = (" FXAT "'{:>7}'.format(at) +
                            "      0      0      0" +
                            '{:>11}'.format('%.4f' % cmd_args.constraintForce) +
                            '{:>11}'.format('%.4f' % cmd_args.constraintWidth) +
                            "     0.0000     0.0000" + "\n")
                    com_file.write(line)
            # BEN
            if len(cmd_args.constrain_by_index):
                line = (" FXTA "'{:>7}'.format(constraint_at_list[0]) +
                                '{:>7}'.format(constraint_at_list[1]) +
                                '{:>7}'.format(constraint_at_list[2]) +
                                '{:>7}'.format(constraint_at_list[3]) + 
                                "  5000.0000   999.0000     0.0000     0.0000\n") 
#  "    99.0000   999.0000     0.0000     0.0000\n")
                com_file.write(line)

        com_file.write(line7)
        com_file.write(line8)
        com_file.write(line9) 
        
    return basename, comfile, inmaefile, outmaefile, logfile


###############################################################################
def get_atom_numbers(coords, st):
    """
    For each coordinate in the list, find the atom that is the closest.

    @type coords:  String of atom coordinates to use.
    @param coords:  str

    @param st:  Current structure.
    @type st:  L{structure.Structure}
    """

    at_list = []
    at_num = None

    coords = coords.split(";")

    if len(coords) != 4:
        print "Could not determine 4 sets of coordinates."
        sys.exit(1)
    else:
        for coord in coords:
            max_dist = 100
            x, y, z = coord.split(",")
            x = float(x.strip(" "))
            y = float(y.strip(" "))
            z = float(z.strip(" "))

            for at in st.atom:
                dist = math.sqrt((x - at.x)**2 + (y - at.y)**2 + (z - at.z)**2)
                if dist < max_dist:
                    max_dist = dist
                    at_num = at.index

            at_list.append(at_num)
    at_list = list(set(at_list))

    if len(at_list) != 4:
        print "Four unique torsion atoms couldn't be identified."
        sys.exit(1)

    return at_list


###############################################################################
def main():
    """
    Main body of the script.
    """

    constraint_at_list = []

    outfile_list = []
    infile_list  = []
    comfile_list = []
    logfile_list = []

    cmd_args = parse_args()
    constraints_enabled = check_args(cmd_args)
  
    if len(cmd_args.constrain_by_index):
         constraint_at_list = cmd_args.constrain_by_index.split(",")

    # BEN     
    basename, comfile, inmaefile, outmaefile, logfile = \
        create_com_file(cmd_args.infile, cmd_args, constraint_at_list)

    outfile_list.append(outmaefile)
    infile_list.append(inmaefile)
    comfile_list.append(comfile)
    logfile_list.append(logfile)
    
    subprocess.call([BMIN, "-NOJOBID", "-LOCAL", basename])
    
    writer = structure.StructureWriter(cmd_args.outfile)

    for out in outfile_list:
        for st in structure.StructureReader(out):
            writer.append(st)
    writer.close()

    file_list = comfile_list + outfile_list + infile_list + logfile_list
   
if __name__ == '__main__':
    cmdline.main_wrapper(main)
