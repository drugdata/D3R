#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_molfilevalidator
----------------------------------

Tests for `molfilevalidator` module.
"""
import sys
import tempfile
import os
import os.path
import shutil
import tarfile

if sys.version_info < (2, 7):
    import unittest2 as unittest
else:
    import unittest

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

from d3r.celpp.task import D3RParameters  # noqa: E402
from d3r import molfilevalidator  # noqa: E402
from d3r.molfilevalidator import D3RAtom  # noqa: E402
from d3r.molfilevalidator import D3RMolecule  # noqa: E402
from d3r.molfilevalidator import D3RMoleculeFromSmileViaOpeneyeFactory \
    # noqa: E402
from d3r.molfilevalidator import D3RMoleculeFromMolFileViaOpeneyeFactory \
    # noqa: E402
from d3r.molfilevalidator import ValidationReport  # noqa: E402
from d3r.molfilevalidator import CompareMolecules  # noqa: E402


class TestMolFileValidator(unittest.TestCase):
    """Tests molfilevalidator.py command line script
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def _get_mol_as_str(self):
        return """
-OEChem-09291710502D

  2  1  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""

    def _get_mol2_as_str(self):
        return """
-OEChem-09291710502D

  2  1  0     0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""

    def test_parse_arguments(self):

        gen_mode = molfilevalidator.GENMOLECULEDB_MODE
        params = molfilevalidator._parse_arguments('hi', [gen_mode])

        self.assertEqual(params.mode, gen_mode)
        self.assertEqual(params.moldir, None)
        self.assertEqual(params.skipligand, None)
        self.assertEqual(params.outputfile, None)
        self.assertEqual(params.moleculedb, None)
        self.assertEqual(params.molcsvligandcol, 0)
        self.assertEqual(params.molcsvsmilecol, 1)
        self.assertEqual(params.loglevel, molfilevalidator.DEFAULT_LOG_LEVEL)
        self.assertEqual(params.excludedir, 'SuppInfo')

        params = molfilevalidator._parse_arguments('hi',
                                                   [gen_mode,
                                                    '--moldir', '/foo',
                                                    '--outputfile', 'out',
                                                    '--skipligand', 'FXR_33',
                                                    '--moleculedb', 'mdb',
                                                    '--log', 'DEBUG',
                                                    '--excludedir', 'hi,how'])
        self.assertEqual(params.mode, gen_mode)
        self.assertEqual(params.moldir, '/foo')
        self.assertEqual(params.skipligand, 'FXR_33')
        self.assertEqual(params.outputfile, 'out')
        self.assertEqual(params.moleculedb, 'mdb')
        self.assertEqual(params.loglevel, 'DEBUG')
        self.assertEqual(params.excludedir, 'hi,how')

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

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
                     'No valid openeye installation or license found')
    def test_d3rmoleculefrommolfileviaopeneyefactory(self):
        mol = self._get_mol_as_str()
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

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
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
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(None)
        self.assertEqual(ha, -1)
        self.assertEqual(mw, -1)
        self.assertEqual(adic, {})
        self.assertEqual(smi_str, '')
        mol = D3RMolecule()

        # test None for atoms
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 0)
        self.assertEqual(mw, 0)
        self.assertEqual(adic, {})
        self.assertEqual(smi_str, '')

        # test with empty list for atoms
        mol.set_atoms([])
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 0)
        self.assertEqual(mw, 0)
        self.assertEqual(adic, {})
        self.assertEqual(smi_str, None)

        # test with 1 atom is hydrogen
        atom = D3RAtom()
        atom.set_is_hydrogen(True)
        atom.set_atomic_number(1)
        mol.set_atoms([atom])
        mol.set_canonical_smiles_str('CC')
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 0)
        self.assertEqual(mw, 0)
        self.assertEqual(adic, {})
        self.assertEqual(smi_str, 'CC')

        # test with 1 atom non hydrogen
        atom.set_is_hydrogen(False)
        atom.set_atomic_number(5)
        mol.set_atoms([atom])
        mol.set_canonical_smiles_str('CCN')
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 1)
        self.assertEqual(mw, 5)
        self.assertEqual(adic, {5: 1})
        self.assertEqual(smi_str, 'CCN')

        # test with 2 atoms where one is 1 non hydrogen
        hatom = D3RAtom()
        hatom.set_is_hydrogen(True)
        hatom.set_atomic_number(1)

        mol.set_atoms([atom, hatom])
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 1)
        self.assertEqual(mw, 5)
        self.assertEqual(adic, {5: 1})
        self.assertEqual(smi_str, 'CCN')

        # test with 4 atoms where one is hydrogen
        catom = D3RAtom()
        catom.set_atomic_number(6)
        catom.set_is_hydrogen(False)
        mol.set_atoms([atom, hatom, catom, catom])
        (ha, mw, adic, smi_str) = molfilevalidator.\
            get_molecule_weight_and_summary(mol)
        self.assertEqual(ha, 3)
        self.assertEqual(mw, 17)
        self.assertEqual(adic, {5: 1, 6: 2})
        self.assertEqual(smi_str, 'CCN')

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
        vr.add_molecule_error('/blah/b.mol', 'XXX_2', (1, 2, {6: 2}, 'CC'),
                              (3, 4, {6: 5}, 'CC'), 'ha')
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
        mol.set_canonical_smiles_str('CC')

        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (3, 5, {5: 1}, 'CC')})
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, False)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 1)
        self.assertEqual(m_errors[0][ValidationReport.MOLFILE], 'm.mol')
        self.assertEqual(m_errors[0][ValidationReport.LIGAND], 'XXX_1')
        self.assertEqual(m_errors[0][ValidationReport.USERMOL],
                         (1, 5, {5: 1}, 'CC'))
        self.assertEqual(m_errors[0][ValidationReport.EXPECTEDMOL],
                         (3, 5, {5: 1}, 'CC'))
        self.assertTrue('Number of heavy atoms '
                        'and or' in m_errors[0][ValidationReport.MESSAGE])

        # test with molecule not matching molecular weight
        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (1, 6, {5: 1}, 'CC')})
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, False)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 1)
        self.assertEqual(m_errors[0][ValidationReport.MOLFILE], 'm.mol')
        self.assertEqual(m_errors[0][ValidationReport.LIGAND], 'XXX_1')
        self.assertEqual(m_errors[0][ValidationReport.USERMOL],
                         (1, 5, {5: 1}, 'CC'))
        self.assertEqual(m_errors[0][ValidationReport.EXPECTEDMOL],
                         (1, 6, {5: 1}, 'CC'))
        self.assertTrue('Number of heavy atoms '
                        'and or' in m_errors[0][ValidationReport.MESSAGE])

        # test with matching molecule but SMILE String differs
        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (1, 5, {5: 1}, 'CN')})
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, False)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 1)
        self.assertEqual(m_errors[0][ValidationReport.MOLFILE], 'm.mol')
        self.assertEqual(m_errors[0][ValidationReport.LIGAND], 'XXX_1')
        self.assertEqual(m_errors[0][ValidationReport.USERMOL],
                         (1, 5, {5: 1}, 'CC'))
        self.assertEqual(m_errors[0][ValidationReport.EXPECTEDMOL],
                         (1, 5, {5: 1}, 'CN'))
        self.assertTrue('Canonical SMILE strings do NOT '
                        'match' in m_errors[0][ValidationReport.MESSAGE])

        # test with matching molecule but SMILE String differs
        # but this time we are skipping smile string comparison
        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (1, 5, {5: 1}, 'CN')},
                              skipsmilecompare=True)
        res = cm.compare_molecules('m.mol', vr, 'XXX_1', mol)
        self.assertEqual(res, True)
        m_errors = vr.get_molecule_errors()
        self.assertEqual(len(m_errors), 0)

        # test with full match
        vr = ValidationReport()
        cm = CompareMolecules({'XXX_1': (1, 5, {5: 1}, 'CC')})
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
            self.assertEqual(str(e), 'Error parsing ligand name from file '
                                     'name: hi')

        res = molfilevalidator._get_ligand_name_from_file_name('DSV-'
                                                               'FXR_12-1.mol')
        self.assertEqual(res, 'FXR_12')

    def test_generate_molecule_database_fromcsv_molcsv_is_none(self):
        params = D3RParameters()
        params.molcsv = None
        res = molfilevalidator.\
            _generate_molecule_database_fromcsv(params, None)
        self.assertEqual(res, 1)

    def test_generate_molecule_database_no_moldir_is_none(self):
        params = D3RParameters()
        params.moldir = None
        res = molfilevalidator.\
            _generate_molecule_database_frommolfiles(params, None)
        self.assertEqual(res, 1)

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
                     'No valid openeye installation or license found')
    def test_generate_molecule_database_from_csv_with_validcsv(self):
        temp_dir = tempfile.mkdtemp()
        try:
            csvfile = os.path.join(temp_dir, 'foo.csv')
            with open(csvfile, 'w') as f:
                f.write('Ligand_ID, Smiles, Target, Subchallenge\n')
                f.write('ABC_1,FC(F),abc,something\n')
                f.write('DEF_1,CC,abc,something\n')
                f.flush()
            params = D3RParameters()
            params.molcsv = csvfile
            params.outputfile = os.path.join(temp_dir, 'res.pickle')
            params.molcsvligandcol = 0
            params.molcsvsmilecol = 1
            factory = D3RMoleculeFromSmileViaOpeneyeFactory()
            res = molfilevalidator._generate_molecule_database_fromcsv(params,
                                                                       factory)
            self.assertEqual(res, 0)
            self.assertTrue(os.path.isfile(params.outputfile))

            params.moleculedb = params.outputfile
            ligand_dic = molfilevalidator._get_molecule_database(params)
            self.assertEqual(len(ligand_dic), 2)
            self.assertEqual(ligand_dic['ABC_1'][0], 3)
            self.assertEqual(ligand_dic['DEF_1'][0], 2)
        finally:
            shutil.rmtree(temp_dir)

    def test_generate_molecule_database_no_molfiles(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.moldir = temp_dir
            output = os.path.join(temp_dir, 'foo.pickle')
            params.outputfile = output
            res = molfilevalidator.\
                _generate_molecule_database_frommolfiles(params, None)
            self.assertEqual(res, 4)
            self.assertFalse(os.path.isfile(output))
        finally:
            shutil.rmtree(temp_dir)

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
                     'No valid openeye installation or license found')
    def test_generate_molecule_database_one_molfile_with_openeye(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.moldir = temp_dir
            output = os.path.join(temp_dir, 'foo.pickle')
            params.outputfile = output
            mymol = os.path.join(temp_dir, 'my-ABC_1-1.mol')
            with open(mymol, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()
            fac = D3RMoleculeFromMolFileViaOpeneyeFactory()
            res = molfilevalidator.\
                _generate_molecule_database_frommolfiles(params, fac)
            self.assertEqual(res, 0)
            self.assertTrue(os.path.isfile(output))

            params.moleculedb = output
            moldb = molfilevalidator._get_molecule_database(params)
            self.assertEqual(len(moldb), 1)
            self.assertEqual(moldb['ABC_1'][0], 2)
            self.assertEqual(moldb['ABC_1'][1], 13)
        finally:
            shutil.rmtree(temp_dir)

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
                     'No valid openeye installation or license found')
    def test_generate_molecule_database_two_molfiles_with_openeye(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.moldir = temp_dir
            output = os.path.join(temp_dir, 'foo.pickle')
            params.outputfile = output
            mymol = os.path.join(temp_dir, 'my-ABC_1-1.mol')
            with open(mymol, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()

            mymol2 = os.path.join(temp_dir, 'yo-XXX_2-3.mol')
            with open(mymol2, 'w') as f:
                f.write(self._get_mol2_as_str())
                f.flush()

            fac = D3RMoleculeFromMolFileViaOpeneyeFactory()
            res = molfilevalidator.\
                _generate_molecule_database_frommolfiles(params, fac)
            self.assertEqual(res, 0)
            self.assertTrue(os.path.isfile(output))

            params.moleculedb = output
            moldb = molfilevalidator._get_molecule_database(params)
            self.assertEqual(len(moldb), 2)
            self.assertEqual(moldb['ABC_1'][0], 2)
            self.assertEqual(moldb['ABC_1'][1], 13)
            self.assertEqual(moldb['XXX_2'][0], 2)
            self.assertEqual(moldb['XXX_2'][1], 14)
        finally:
            shutil.rmtree(temp_dir)

    def test_molfile_from_tarfile_generator_no_tarfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            n = os.path.join(temp_dir, 'foo.tar.gz')
            try:
                list(molfilevalidator._molfile_from_tarfile_generator(n))
                self.fail('Expected IOError')
            except IOError as e:
                self.assertTrue('No such file' in str(e))

        finally:
            shutil.rmtree(temp_dir)

    def test_molfile_from_tarfile_generator_no_mols_in_tar(self):
        temp_dir = tempfile.mkdtemp()
        try:
            afile = os.path.join(temp_dir, 'foo.pdb')
            with open(afile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()
            myf = os.path.join(temp_dir, 'foo.tar.gz')
            tf = tarfile.open(myf, 'w:gz')
            tf.add(afile)
            tf.close()

            mols = list(molfilevalidator._molfile_from_tarfile_generator(myf))
            self.assertEqual(len(mols), 0)

            mols = list(molfilevalidator.
                        _molfile_from_tarfile_generator(myf,
                                                        direxclude='blah'))
            self.assertEqual(len(mols), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_molfile_from_tarfile_generator_one_mol_in_tar(self):
        temp_dir = tempfile.mkdtemp()
        try:

            afile = os.path.join(temp_dir, 'foo.pdb')
            with open(afile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()

            mfile = os.path.join(temp_dir, 'yo.mol')
            with open(mfile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()

            myf = os.path.join(temp_dir, 'foo.tar.gz')
            tf = tarfile.open(myf, 'w:gz')
            tf.add(afile)
            tf.add(mfile, arcname='SuppInfo/hi.mol')
            tf.close()

            mols = list(molfilevalidator._molfile_from_tarfile_generator(myf))
            self.assertEqual(len(mols), 1)
            self.assertTrue(mols[0].endswith('hi.mol'))

            # try with exclude enabled not match
            mols = list(molfilevalidator.
                        _molfile_from_tarfile_generator(myf,
                                                        direxclude='Supp,yo'))
            self.assertEqual(len(mols), 1)

            # try with exclude enabled
            mols = list(molfilevalidator.
                        _molfile_from_tarfile_generator(myf,
                                                        direxclude='SuppInfo'))
            self.assertEqual(len(mols), 0)
        finally:
            shutil.rmtree(temp_dir)

    def test_molfile_from_tarfile_generator_two_mol_in_tar(self):
        temp_dir = tempfile.mkdtemp()
        try:
            afile = os.path.join(temp_dir, 'foo.pdb')
            with open(afile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()

            mfile = os.path.join(temp_dir, 'yo.mol')
            with open(mfile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()

            m2file = os.path.join(temp_dir, 'xxd.mol')
            with open(m2file, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()

            myf = os.path.join(temp_dir, 'foo.tar.gz')
            tf = tarfile.open(myf, 'w:gz')
            tf.add(afile)
            tf.add(mfile, arcname='hi.mol')
            tf.add(m2file, arcname='bye.mol')
            tf.close()

            mols = list(molfilevalidator._molfile_from_tarfile_generator(myf))
            self.assertEqual(len(mols), 2)
            self.assertTrue(mols[0].endswith('hi.mol'))
            self.assertTrue(mols[1].endswith('bye.mol'))
        finally:
            shutil.rmtree(temp_dir)

    def test_validate_molfiles_in_tarball_usersubmission_is_none(self):
        params = D3RParameters()
        params.usersubmission = None
        res = molfilevalidator._validate_molfiles_in_tarball(params, None,
                                                             None)
        self.assertEqual(res, 1)

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
                     'No valid openeye installation or license found')
    def test_validate_molfiles_in_tarball_real(self):
        temp_dir = tempfile.mkdtemp()
        try:

            mfile = os.path.join(temp_dir, 'yo.mol')
            with open(mfile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()
            myf = os.path.join(temp_dir, 'foo.tar.gz')

            bfile = os.path.join(temp_dir, 'x.mol')
            with open(bfile, 'w') as f:
                f.write('blah\n')
                f.flush()
                f.close()
            tf = tarfile.open(myf, 'w:gz')
            tf.add(mfile, arcname='yo.mol')
            tf.add(mfile, arcname='x-ABC_1-1.mol')
            tf.add(bfile, arcname='z-BBB_2-213.mol')
            tf.close()

            params = D3RParameters()
            params.usersubmission = myf
            params.excludedir = 'SuppInfo'
            moldb = {}
            moldb['ABC_1'] = (1, 2, {7: 2}, 'CC')
            molfactory = D3RMoleculeFromMolFileViaOpeneyeFactory()
            res = molfilevalidator._validate_molfiles_in_tarball(params,
                                                                 molfactory,
                                                                 moldb)
            self.assertTrue('parsing ligand name from file '
                            'name: yo.mol' in res.get_as_string())
            self.assertTrue('1 non hydrogen atoms, but '
                            'got 2' in res.get_as_string())
            self.assertTrue('2 for non hydrogen atomic weight, but '
                            'got 13' in res.get_as_string())
            self.assertTrue('BBB_2\n\tligand not in molecule '
                            in res.get_as_string())
        finally:
            shutil.rmtree(temp_dir)

    @unittest.skipIf(OPENEYE_MOD_LOADED is False,
                     'No valid openeye installation or license found')
    def test_validate_molfiles_in_tarball_real_with_skiplist(self):
        temp_dir = tempfile.mkdtemp()
        try:
            myf = os.path.join(temp_dir, 'foo.tar.gz')
            mfile = os.path.join(temp_dir, 'yo.mol')
            with open(mfile, 'w') as f:
                f.write(self._get_mol_as_str())
                f.flush()
            with tarfile.open(myf, 'w:gz') as tf:

                tf.add(mfile, arcname='x-ABC_1-1.mol')

            params = D3RParameters()
            params.usersubmission = myf
            params.excludedir = 'SuppInfo'
            params.skipligand = 'XXX_3,ABC_1'
            moldb = {}
            moldb['ABC_1'] = (1, 2, {7: 2})
            molfactory = D3RMoleculeFromMolFileViaOpeneyeFactory()
            res = molfilevalidator._validate_molfiles_in_tarball(params,
                                                                 molfactory,
                                                                 moldb)
            self.assertEqual(res.get_as_string(), '')
        finally:
            shutil.rmtree(temp_dir)

    def test_main_with_error(self):
        res = molfilevalidator.main(['yo', molfilevalidator.VALIDATE_MODE])
        self.assertEqual(res, 2)

        res = molfilevalidator.main(['yo',
                                     molfilevalidator.GENMOLECULEDB_MODE])
        self.assertEqual(res, 3)


if __name__ == '__main__':
    unittest.main()
