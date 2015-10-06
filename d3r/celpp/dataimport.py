# -*- coding: utf-8 -*-

__author__ = 'churas'

import logging
import os
from d3r.celpp.task import D3RTask

logger = logging.getLogger(__name__)


class DataImportTask(D3RTask):
    """Represents DataImport Task

       """
    NONPOLYMER_TSV = "new_release_structure_nonpolymer.tsv"
    SEQUENCE_TSV = "new_release_structure_sequence.tsv"
    CRYSTALPH_TSV = "new_release_crystallization_pH.tsv"

    def __init__(self, path, args):
        super(DataImportTask, self).__init__(path, args)
        self.set_name('dataimport')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def get_nonpolymer_tsv(self):
        """Returns path to new_release_structure_nonpolymer.tsv file
        :return: full path to DataImportTask.NONPOLYMER_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.NONPOLYMER_TSV)

    def get_sequence_tsv(self):
        """Returns path to new_release_structure_sequence.tsv file
        :return: full path to DataImportTask.SEQUENCE_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.SEQUENCE_TSV)

    def get_crystalph_tsv(self):
        """Returns path to new_release_crystallization_pH.tsv file
        :return: full path to DataImportTask.CRYSTALPH_TSV file
        """
        return os.path.join(self.get_dir(),
                            DataImportTask.CRYSTALPH_TSV)
