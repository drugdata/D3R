# -*- coding: utf-8 -*-

import os
import logging
import re
logger = logging.getLogger(__name__)


def find_latest_year(celppdir):
    """Given a directory find the latest year

       The latest year will be a folder of 4 digits
       and have the highest value ie 2015
       :return: Directory or None if none is found
       :raises: Exception: If celppdir is not a directory
       """

    if not os.path.isdir(celppdir):
        raise Exception(celppdir+" is not a directory")

    dir_pattern = re.compile("^\d\d\d\d$")

    latest_year = -1
    latest_entry = None
    full_path = None
    for entry in os.listdir(celppdir):
        if re.match(dir_pattern, entry):
            entry_year = int(entry)

            if entry_year > latest_year:
                full_path = os.path.join(celppdir, entry)
                if os.path.isdir(full_path):
                    latest_year = entry_year
                    latest_entry = full_path

    if latest_year == -1:
        return

    return latest_entry


def find_latest_weekly_dataset(celppdir):
    """Given a directory find the latest dataset

       This method looks in directory passed in for
       paths with dataset.week.# format and returns
       the path with highest # on it.  If none are
       found then None is returned.
       :return: Directory upon success or None if no directory found
       :raises: Exception: If celppdir is not a directory
       """

    latest_year = find_latest_year(celppdir)

    if latest_year is None:
        return

    dir_pattern = re.compile("^dataset.week.\d+$")

    latest_entry = None
    latest_weekno = -1
    full_path = None
    for entry in os.listdir(latest_year):
        if re.match(dir_pattern, entry):
            weekno = re.sub("dataset.week.", "", entry)
            if weekno > latest_weekno:
                full_path = os.path.join(latest_year, entry)
                if os.path.isdir(full_path):
                    latest_weekno = weekno
                    latest_entry = full_path

    return latest_entry
