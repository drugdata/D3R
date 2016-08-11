__author__ = 'churas'

"""
test_vinadocking
--------------------------------

Tests for `vinadocking` module.
"""

import unittest
import os
import os.path
import re


class TestCelppCommandLineScripts(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def get_list_of_celpp_scripts(self):
        d3r_dir = os.path.join(os.getcwd(), 'd3r')
        celppscripts = []
        for entry in os.listdir(d3r_dir):
            fullpath = os.path.join(d3r_dir, entry)
            if fullpath.endswith('__init__.py'):
                continue
            if os.path.isfile(fullpath):
                if fullpath.endswith('.py'):
                    celppscripts.append(fullpath)

        return celppscripts

    def test_every_script_starts_with_correct_shebang(self):
        """This test was added since a commit was recently made
           where the shebang was modified and had a typo in it
           causing all the scripts to fail
        """
        for script in self.get_list_of_celpp_scripts():
            f = open(script, 'r')
            val = re.sub('\n$', '', f.readline())
            self.assertTrue(re.match('#! *\/usr\/bin\/env +python', val),
                            ' invalid shebang line in file: ' + script +
                            ' : ' + val)
            f.close()


if __name__ == '__main__':
    unittest.main()
