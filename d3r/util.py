# -*- coding: utf-8 -*-

import os
import logging
import re
from datetime import date
from datetime import timedelta

logger = logging.getLogger(__name__)

DAYS_IN_WEEK = 7


def find_latest_year(celppdir):
    """Given a directory find the latest year

       The latest year will be a folder of 4 digits
       and have the highest value ie 2015
       :return: Directory or None if none is found
       :raises: Exception: If celppdir is not a directory
       """

    if not os.path.isdir(celppdir):
        raise Exception(celppdir + " is not a directory")

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


def get_celpp_week_of_year_from_date(the_date):
    """Returns CELPP week of year
    The CELPP week of year matches ISO week of year, but it
    starts on Friday instead of Monday.  The way it works is
    the code examines `the_date` and if its 12:00am Monday -
    11:59pm Thursday the current week according to date.isocalendar()
    is used.  If date falls between Friday 12:00am to Sunday 11:59pm
    this function advances the date to the following Monday and extracts
    week of year.
       :return: tuple with first value the week and second value the year
       :raises: Exception if `the_date` is None
    """
    if the_date is None:
        raise Exception('Must pass a valid date')
    assert isinstance(the_date, date)
    updated_date = the_date

    day_of_week = the_date.weekday()
    if day_of_week >= 4:
        updated_date = the_date + timedelta(days=(DAYS_IN_WEEK - day_of_week))
        logger.debug('Date ' + str(the_date) + ' is between Friday and ' +
                     'Sunday shifting date to ' + str(updated_date))

    return (updated_date.isocalendar()[1],
            updated_date.isocalendar()[0])

def create_celpp_week_dir(celpp_week_tuple, celppdir):
    """Creates celpp week dir in `celppdir`
    Creates celpp week directory under `celppdir`
    using `celpp_week_tuple` to get year and week of year
    diretories created:
    celppdir/<celpp_week_tuple[1]>/dataset.week.<celpp_week_tuple[0]>
    If directory already exists no operation is performed.
    :raises: OSError if unable to create directory
    """
    if celpp_week_tuple is None:
        raise Exception('celpp week tuple must be set')

    if celppdir is None:
        raise Exception('celppdir must be set')

    dir_to_create = os.path.join(celppdir, str(celpp_week_tuple[1]),
                                 'dataset.week.'+str(celpp_week_tuple[0]))
    if os.path.isdir(dir_to_create):
        logger.debug(dir_to_create + ' directory already exists skipping creation')
        return

    logger.debug('Creating ' + dir_to_create)
    os.makedirs(dir_to_create, 0775)
