# -*- coding: utf-8 -*-

__author__ = 'churas'

import os
import logging
from d3r.celpp.task import D3RTask

logger = logging.getLogger(__name__)


class MakeBlastDBTask(D3RTask):
    """Represents Blastable database

       This object does not follow the standard D3RTask
       structure.  This is because the blastable database
       is created externally in a legacy structure.

       """

    def __init__(self, path, args):
        """Constructor

           Constructor sets name to makeblastdb and ignores
           path set in <path> variable and instead sets
           path to args.blastdir. stage is set to 1

        """
        self.set_args(args)
        self.set_path(args.blastdir)
        self.set_name('makeblastdb')
        self.set_stage(1)
        self.set_status(D3RTask.UNKNOWN_STATUS)

    def create_dir(self):
        """Creates path set in get_path()

           """
        the_path = os.path.join(self._path, self.get_dir_name())
        os.makedirs(the_path)
        return the_path

    def get_dir_name(self):
        """Will always return current
        :return: 'current' string
        """
        return 'current'
