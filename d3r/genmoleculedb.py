#! /usr/bin/env python

import sys
import os
import argparse
import glob
import logging
from openeye import oechem
import pickle

import d3r
from d3r.celpp import util
from d3r.celpp.task import D3RParameters


# create logger
logger = logging.getLogger('d3r.genmoleculedb')
DEFAULT_LOG_LEVEL = 'ERROR'
p = D3RParameters()
p.loglevel = DEFAULT_LOG_LEVEL
util.setup_logging(p)


def _parse_arguments(desc, args):
    """Parses command line arguments using argparse.
    """
    parsed_arguments = D3RParameters()

    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("moldir", help='Directory containing .mol files')
    parser.add_argument("outputdb", help='Output database file')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' + d3r.__version__))
    return parser.parse_args(args, namespace=parsed_arguments)


def get_mw_na (rdmol):
    #from rdmol object get the molecular weight and number of heavy atoms
    atom_dic = {}
    heavy_atom = 0
    molecular_weight = 0
    for atom in rdmol.GetAtoms():
        if not atom.IsHydrogen():
            heavy_atom +=1
            atomical_number = atom.GetAtomicNum()
            molecular_weight += atomical_number
            logger.debug('Atom name: ' + atom.GetName())
            if atomical_number not in atom_dic:
                atom_dic[atomical_number] = 1
            else:
                atom_dic[atomical_number] += 1
    return heavy_atom, molecular_weight, atom_dic


def _generate_molecule_database(theargs):
    search_path = os.path.join(theargs.moldir, '*.mol')
    all_mol_files = glob.glob(search_path)
    ligand_dic = {}
    for mol_file in all_mol_files:
        logger.debug('File: ' + mol_file)
        ligand_name = mol_file.split("-")[1]
        logger.debug('Ligand: ' + ligand_name)
        istream = oechem.oemolistream()
        istream.open(mol_file)
        openeye_mol = oechem.OEMol()
        oechem.OEReadMolecule(istream, openeye_mol)
        istream.close()
        if not ligand_name in ligand_dic:
            ligand_dic[ligand_name] = get_mw_na(openeye_mol)

    logger.debug('Writing output to: ' + theargs.outputdb)
    p_f = open(theargs.outputdb, "w")
    pickle.dump(ligand_dic, p_f)
    p_f.close()
    return 0


def main(args):
    """Main entry into genmoleculedb
    :param args: should be set to sys.argv which is a list of arguments
                 starting with script name as the first argument
    """
    desc = """
              Version {version}

              Parses a set of .mol files found in input directory and
              generates a database used by molfilevalidator.py
              (http://www.drugdesigndata.org)

              """.format(version=d3r.__version__)

    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = d3r.__version__

    util.setup_logging(theargs)

    try:
        return _generate_molecule_database(theargs)
    except Exception:
        logger.exception("Error caught exception")
        return 2


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
