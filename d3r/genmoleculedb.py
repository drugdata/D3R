#! /usr/bin/env python

import sys
import os
import argparse
import glob
import logging
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

LIG_NAME_DELIM = '-'


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
    parser.add_argument("--log", dest="loglevel",
                        choices=['DEBUG', 'INFO', DEFAULT_LOG_LEVEL,
                                 'ERROR', 'CRITICAL'],
                        help="Set the logging level (default " +
                        DEFAULT_LOG_LEVEL + ")", default=DEFAULT_LOG_LEVEL)
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
            logger.debug('Atom (' + atom.GetName() + ') atomic #: ' +
                         str(atomical_number))
            if atomical_number not in atom_dic:
                atom_dic[atomical_number] = 1
            else:
                atom_dic[atomical_number] += 1
    return heavy_atom, molecular_weight, atom_dic


def _get_ligand_name_from_file_name(file_name):
    """Extracts ligand name from file name. It is assume
    the file name has format of XXXX-<LIGAND NAME>-X.mol
    """
    if file_name is None:
        raise ValueError('file_name cannot be None')

    base_fname = os.path.basename(file_name)
    hyphen_split = base_fname.split(LIG_NAME_DELIM)
    if len(hyphen_split) < 2:
        raise ValueError('Error parsing ligand name from file name: ' +
                         base_fname)
    return hyphen_split[1]


def _generate_molecule_database(theargs):
    search_path = os.path.join(theargs.moldir, '*.mol')
    all_mol_files = glob.glob(search_path)
    ligand_dic = {}
    for mol_file in all_mol_files:
        logger.info('Reading file: ' + mol_file)
        ligand_name = _get_ligand_name_from_file_name(mol_file)
        logger.info('Ligand: ' + ligand_name)

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
        from openeye import oechem
        return _generate_molecule_database(theargs)
    except Exception:
        logger.exception("Error caught exception")
        return 2


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
