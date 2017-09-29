#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_molfilevalidator
----------------------------------

Tests for `molfilevalidator` module.
"""

import unittest
import tempfile
import os
import os.path
import sys
import shutil

try:
    OPENEYE_MOD_LOADED = False
    from openeye import oechem
    try:
        if oechem.OEChemIsLicensed():
            OPENEYE_MOD_LOADED = True
        else:
            sys.stderr.write('WARNING: No valid license found for openeye ' +
                             'Please verify OE_LICENSE environment variable'
                             'is set to valid license file\n')
    except AttributeError as ae:
        sys.stderr.write('WARNING Unable to check if openey is licensed' +
                         str(ae))
except ImportError as e:
    pass

from d3r import molfilevalidator
from d3r.molfilevalidator import D3RAtom
from d3r.molfilevalidator import D3RMolecule
from d3r.molfilevalidator import D3RMoleculeFromSmileViaOpeneyeFactory
from d3r.molfilevalidator import D3RMoleculeFromMolFileViaOpeneyeFactory
from d3r.molfilevalidator import ValidationReport
from d3r.molfilevalidator import CompareMolecules


class TestMolFileValidator(unittest.TestCase):
    """Tests molfilevalidator.py command line script
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_parse_arguments(self):

        gen_mode = molfilevalidator.GENMOLECULEDB_MODE
        params = molfilevalidator._parse_arguments('hi', [gen_mode])

        self.assertEqual(params.mode, gen_mode)
        self.assertEqual(params.moldir, None)
        self.assertEqual(params.skipligand, None)
        self.assertEqual(params.outputfile, None)
        self.assertEqual(params.moleculedb, None)
        self.assertEqual(params.loglevel, molfilevalidator.DEFAULT_LOG_LEVEL)

        params = molfilevalidator._parse_arguments('hi',
                                                   [gen_mode,
                                                    '--moldir', '/foo',
                                                    '--outputfile', 'out',
                                                    '--skipligand', 'FXR_33',
                                                    '--moleculedb', 'mdb',
                                                    '--log', 'DEBUG'])
        self.assertEqual(params.mode, gen_mode)
        self.assertEqual(params.moldir, '/foo')
        self.assertEqual(params.skipligand, 'FXR_33')
        self.assertEqual(params.outputfile, 'out')
        self.assertEqual(params.moleculedb, 'mdb')
        self.assertEqual(params.loglevel, 'DEBUG')

        val_mode = molfilevalidator.VALIDATE_MODE
        params = molfilevalidator._parse_arguments('hi', [val_mode])
        self.assertEqual(params.mode, val_mode)

    def test_d3ratom(self):
        atom = D3RAtom()
        self.assertEqual(atom.is_hydrogen(), False)
        self.assertEqual(atom.get_atomic_number(), 0)
        self.assertEqual(atom.get_atomic_name(), '')

        atom.set_is_hydrogen(True)
        atom.set_atomic_number(3)
        atom.set_atomic_name('Lithium')

        self.assertEqual(atom.is_hydrogen(), True)
        self.assertEqual(atom.get_atomic_number(), 3)
        self.assertEqual(atom.get_atomic_name(), 'Lithium')

        atom.set_is_hydrogen(False)
        self.assertEqual(atom.is_hydrogen(), False)

    def test_d3rmolecule(self):
        mol = D3RMolecule()
        self.assertEqual(mol.get_atoms(), None)
        mol.set_atoms(['hi', 'bye'])
        self.assertEqual(mol.get_atoms(), ['hi', 'bye'])

    @unittest.skipIf(OPENEYE_MOD_LOADED == False,
                     'No valid openeye installation or license found')
    def test_d3rmoleculefrommolfileviaopeneyefactory(self):

        mol = """
-OEChem-09291710502D

  2  1  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""
        temp_dir = tempfile.mkdtemp()
        try:
            factory = D3RMoleculeFromMolFileViaOpeneyeFactory()
            badmol = os.path.join(temp_dir, 'bad.mol')

            with open(badmol, 'w') as f:
                f.write('hi\n')
                f.flush()
                f.close()
            try:
                factory.get_d3rmolecule(badmol)
                self.fail('Expected ValueError')
            except ValueError as e:
                self.assertEqual(str(e), 'OEReadMolecule returned False when '
                                         'trying to read mol file')
            mfile = os.path.join(temp_dir, 'foo.mol')
            with open(mfile, 'w') as f:
                f.write(mol)
                f.flush()
                f.close()
            d3rmol = factory.get_d3rmolecule(mfile)
            self.assertEqual(len(d3rmol.get_atoms()), 2)
            self.assertEqual(d3rmol.get_atoms()[0].is_hydrogen(), False)
            self.assertEqual(d3rmol.get_atoms()[0].get_atomic_number(), 6)
            self.assertEqual(d3rmol.get_atoms()[1].is_hydrogen(), False)
            self.assertEqual(d3rmol.get_atoms()[1].get_atomic_number(), 7)
        finally:
            shutil.rmtree(temp_dir)

    @unittest.skipIf(OPENEYE_MOD_LOADED == False,
                     'No valid openeye installation or license found')
    def test_d3rmoleculefromsmileviaopeneyefactory(self):
        factory = D3RMoleculeFromSmileViaOpeneyeFactory()
        d3rmol = factory.get_d3rmolecule('CN')
        self.assertEqual(len(d3rmol.get_atoms()), 2)
        self.assertEqual(d3rmol.get_atoms()[0].is_hydrogen(), False)
        self.assertEqual(d3rmol.get_atoms()[0].get_atomic_number(), 6)
        self.assertEqual(d3rmol.get_atoms()[1].is_hydrogen(), False)
        self.assertEqual(d3rmol.get_atoms()[1].get_atomic_number(), 7)


    def test_get_molecule_weight_and_summary(self):
        # test where molecule is None
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(None)
        self.assertEqual(ha, -1)
        self.assertEqual(mw, -1)
        self.assertEqual(adic, {})
        mol = D3RMolecule()

        # test None for atoms
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 0)
        self.assertEqual(mw, 0)
        self.assertEqual(adic, {})

        # test with empty list for atoms
        mol.set_atoms([])
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 0)
        self.assertEqual(mw, 0)
        self.assertEqual(adic, {})

        # test with 1 atom is hydrogen
        atom = D3RAtom()
        atom.set_is_hydrogen(True)
        atom.set_atomic_number(1)
        mol.set_atoms([atom])
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 0)
        self.assertEqual(mw, 0)
        self.assertEqual(adic, {})

        # test with 1 atom non hydrogen
        atom.set_is_hydrogen(False)
        atom.set_atomic_number(5)
        mol.set_atoms([atom])
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 1)
        self.assertEqual(mw, 5)
        self.assertEqual(adic, {5: 1})

        # test with 2 atoms where one is 1 non hydrogen
        hatom = D3RAtom()
        hatom.set_is_hydrogen(True)
        hatom.set_atomic_number(1)

        mol.set_atoms([atom, hatom])
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 1)
        self.assertEqual(mw, 5)
        self.assertEqual(adic, {5: 1})

        # test with 4 atoms where one is hydrogen
        catom = D3RAtom()
        catom.set_atomic_number(6)
        catom.set_is_hydrogen(False)
        mol.set_atoms([atom, hatom, catom, catom])
        (ha, mw, adic) = molfilevalidator.get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 3)
        self.assertEqual(mw, 17)
        self.assertEqual(adic, {5: 1, 6: 2})

    def test_validation_report(self):
        vr = ValidationReport()

        # test on empty
        self.assertEqual(vr.get_as_string(), '')

        # add a ligand error with none for all fields
        vr.add_ligand_error(None, None, None)
        self.assertTrue('\nLigand Errors\n'
                        '------------\n\n' in vr.get_as_string())
        self.assertTrue('\n\nIn file: None ligand: '
                        'None\n\tNone' in vr.get_as_string())

        # add a molecule with no tuples
        vr = ValidationReport()
        vr.add_molecule_error('x.mol', 'DDD_3', None, None, 'yo')
        self.assertTrue('\nMolecule Errors\n' in vr.get_as_string())
        self.assertTrue('\n\nIn file: x.mol ligand: '
                        'DDD_3\n\tyo\n\n' in vr.get_as_string())

        vr = ValidationReport()
        vr.add_ligand_error('/foo/my.mol', 'ABC_1', 'some error')
        self.assertTrue('\nLigand Errors\n'
                        '------------\n\n' in vr.get_as_string())
        self.assertTrue('\n\nIn file: my.mol ligand: '
                        'ABC_1\n\tsome error' in vr.get_as_string())

        # add a molecule error
        vr.add_molecule_error('/blah/b.mol', 'XXX_2', (1, 2, {6: 2}),
                              (3, 4, {6: 5}), 'ha')
        self.assertTrue('\nLigand Errors\n'
                        '------------\n\n' in vr.get_as_string())
        self.assertTrue('\n\nIn file: my.mol ligand: '
                        'ABC_1\n\tsome error' in vr.get_as_string())
        self.assertTrue('\nMolecule Errors\n' in vr.get_as_string())

        self.assertTrue('\n\nIn file: b.mol ligand: '
                        'XXX_2 ha' in vr.get_as_string())
        self.assertTrue('\n\tExpected 3 non hydrogen atoms, '
                         'but got 1' in vr.get_as_string())
        self.assertTrue('\n\tExpected 4 for non hydrogen atomic weight, '
                         'but got 2' in vr.get_as_string())
        self.assertTrue('\n\tExpected atom map { atomic #: # atoms,...} '
                        '{6: 5}, but got {6: 2}\n\n' in vr.get_as_string())


    def test_compare_molecules(self):
        cm = CompareMolecules({})
        vr = ValidationReport()

        # test with molecule not in db
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', D3RMolecule())
        self.assertEqual(res, False)
        self.assertEqual(len(vr.get_molecule_errors()), 0)
        ligand_errors = vr.get_ligand_errors()
        self.assertEqual(len(ligand_errors), 1)
        self.assertEqual(ligand_errors[0][ValidationReport.MOLFILE], 'm.mol')
        self.assertEqual(ligand_errors[0][ValidationReport.LIGAND], 'XXX_1')
        self.assertEqual(ligand_errors[0][ValidationReport.MESSAGE],
                         'ligand not in molecule database')

        # test with molecule not matching number heavy atoms
        mol = D3RMolecule()
        atom = D3RAtom()
        atom.set_is_hydrogen(False)
        atom.set_atomic_number(5)
        mol.set_atoms([atom])

        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (3, 5, {5: 1})})
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, False)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 1)
        self.assertEqual(m_errors[0][ValidationReport.MOLFILE], 'm.mol')
        self.assertEqual(m_errors[0][ValidationReport.LIGAND], 'XXX_1')
        self.assertEqual(m_errors[0][ValidationReport.USERMOL], (1, 5, {5: 1}))
        self.assertEqual(m_errors[0][ValidationReport.EXPECTEDMOL],
                         (3, 5, {5: 1}))
        self.assertTrue('Number of heavy '
                        'atoms and or' in m_errors[0][ValidationReport.MESSAGE])

        # test with molecule not matching molecular weight
        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (1, 6, {5: 1})})
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, False)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 1)
        self.assertEqual(m_errors[0][ValidationReport.MOLFILE], 'm.mol')
        self.assertEqual(m_errors[0][ValidationReport.LIGAND], 'XXX_1')
        self.assertEqual(m_errors[0][ValidationReport.USERMOL], (1, 5, {5: 1}))
        self.assertEqual(m_errors[0][ValidationReport.EXPECTEDMOL],
                         (1, 6, {5: 1}))
        self.assertTrue('Number of heavy '
                        'atoms and or' in m_errors[0][ValidationReport.MESSAGE])

        # test with matching molecule
        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (1, 5, {5: 1})})
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, True)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 0)

    def test_get_ligand_name_from_file_name(self):
        try:
            molfilevalidator._get_ligand_name_from_file_name(None)
            self.fail('Expected ValueError')
        except ValueError as e:
            self.assertEqual(str(e), 'file_name cannot be None')

        try:
            molfilevalidator._get_ligand_name_from_file_name('hi')
            self.fail('Expected ValueError')
        except ValueError as e:
            self.assertEqual(str(e), 'Error parsing ligand name from file name: hi')

        res = molfilevalidator._get_ligand_name_from_file_name('DSV-FXR_12-1.mol')
        self.assertEqual(res, 'FXR_12')



if __name__ == '__main__':
    unittest.main()
