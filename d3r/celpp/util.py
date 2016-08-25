# -*- coding: utf-8 -*-

import os
import logging
import re
import urllib
import gzip
import subprocess
import shlex
import time
import urllib2
from dateutil.parser import parse
from datetime import date
from datetime import datetime
from datetime import timedelta
from dateutil.tz import tzlocal


logger = logging.getLogger(__name__)

DAYS_IN_WEEK = 7

FRIDAY_WEEKDAY = 4

DATA_SET_WEEK_PREFIX = 'dataset.week.'

LOG_FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"


class DownloadError(Exception):
    """Exception to denote error downloading data
    """
    pass


def get_all_celpp_years(celppdir):
    """Given a directory find all year directories within

       This method looks for all directories within `celppdir`
       that fit the structure of a year directory ####
       or the regex ^\d\d\d\d$.
       Full paths to those directories are returned
       :return: List year directories, note path is NOT prefixed
       :raises: Exception if celppdir is not a directory
    """
    if not os.path.isdir(celppdir):
        raise Exception(celppdir + " is not a directory")

    year_list = []
    dir_pattern = re.compile("^\d\d\d\d$")
    for entry in os.listdir(celppdir):
        if re.match(dir_pattern, entry):
            full_path = os.path.join(celppdir, entry)
            if os.path.isdir(full_path):
                year_list.append(entry)

    return year_list


def find_latest_year(celppdir):
    """Given a directory find the latest year

       The latest year will be a folder of 4 digits
       and have the highest value ie 2015
       :return: Directory or None if none is found
       :raises: Exception: If celppdir is not a directory
       """

    if not os.path.isdir(celppdir):
        raise Exception(celppdir + " is not a directory")

    latest_year = -1
    latest_entry = None
    full_path = None
    for entry in get_all_celpp_years(celppdir):
            entry_year = int(entry)
            if entry_year > latest_year:
                full_path = os.path.join(celppdir, entry)
                if os.path.isdir(full_path):
                    latest_year = entry_year
                    latest_entry = full_path

    if latest_year == -1:
        return

    return latest_entry


def get_celpp_year_from_path(dir):
    """Given path extract the year number from it

       Expects path of
       /..../year/dataset.week.#
       or year/dataset.week.#

    :param dir: full or partial path to week
    :return year #### from path or '0' if year extracted is not a 4
    digit number
    """
    if dir is None:
        raise Exception('Path is none')

    remove_week = re.sub("\/" + DATA_SET_WEEK_PREFIX + ".*$", "", dir)
    year = re.sub("^.*\/", "", remove_week)

    dir_pattern = re.compile("^\d\d\d\d$")

    if not re.match(dir_pattern, year):
        logger.warning('Year ( ' + year + ' ) extracted from ' + dir +
                       ' is not a 4 digit number')
        return '0'
    return year


def get_all_celpp_weeks(celpp_year_dir):
    """Given a celpp year directory get all week directories

       Expects celpp_year_dir to be full path to year dir ie ####/
       and finds all directories of format dataset.week.#
       returning a list of the week numbers ie #
       :returns: List of week numbers
       :raises Exception: If celpp_year_dir is not a directory
    """
    if not os.path.isdir(celpp_year_dir):
        raise Exception(celpp_year_dir + " is not a directory")

    week_list = []

    dir_pattern = re.compile("^" + DATA_SET_WEEK_PREFIX + "\d+$")

    for entry in os.listdir(celpp_year_dir):
        if re.match(dir_pattern, entry):
            weeknum = re.sub(DATA_SET_WEEK_PREFIX, "", entry)
            full_path = os.path.join(celpp_year_dir, entry)
            if os.path.isdir(full_path):
                week_list.append(weeknum)

    return week_list


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

    latest_entry = None
    latest_weeknum = -1
    full_path = None
    for entry in get_all_celpp_weeks(latest_year):
            weeknum = int(entry)
            if weeknum > latest_weeknum:
                full_path = os.path.join(latest_year,
                                         DATA_SET_WEEK_PREFIX + entry)
                if os.path.isdir(full_path):
                    latest_weeknum = weeknum
                    latest_entry = full_path

    return latest_entry


def get_celpp_week_number_from_path(dir):
    """Given path extract the week number from it

       Expects either path of
       /..../year/dataset.week.#
       or dataset.week.#

    :param dir: full or partial path to week
    :return the # from dataset.week.# in `dir` path or '0' if value is
    not a 1-2 digit number
    """
    if dir is None:
        raise Exception('Path is none')

    weeknum = re.sub("^.*" + DATA_SET_WEEK_PREFIX, "", dir)

    dir_pattern = re.compile("^\d\d?$")
    if not re.match(dir_pattern, weeknum):
        logger.warning('Week number ( ' + weeknum + ' ) extracted from path '
                       + dir + ' is not a 2 digit number')
        return '0'
    return weeknum


def get_previous_friday_from_date(the_datetime):
    """Given a date, function finds date of previous Friday
       Given `the_date` code finds the date of the previous
       Friday with same hour, minute, and second.
       :return: datetime object containing date of previous Friday
       :raises: Exception if `the_date` is None
    """
    if the_datetime is None:
        raise Exception('Must pass a valid datetime')
    assert isinstance(the_datetime, datetime)

    weekday = the_datetime.weekday()
    if weekday == FRIDAY_WEEKDAY:
        return the_datetime

    if weekday > FRIDAY_WEEKDAY:
        tdelta = timedelta(days=(weekday - FRIDAY_WEEKDAY))
    else:
        if weekday < FRIDAY_WEEKDAY:
            tdelta = timedelta(days=(3 + weekday))

    return the_datetime - tdelta


def is_datetime_after_celpp_week_start(the_datetime):
    """Determines if `the_datetime` is after celpp week start
    :returns: True if yes otherwise False
    :raises Exception: If `the_datetime` is not a valid datetime object
    """
    if the_datetime is None:
        raise Exception('Must pass a valid datetime')
    assert isinstance(the_datetime, datetime)

    prev_friday = get_previous_friday_from_date(datetime.now(tzlocal()))

    if the_datetime >= prev_friday:
        return True

    return False


def has_url_been_updated_since_start_of_celpp_week(url):
    """Checks if url has been updated since start of celpp week
    """
    if url is None:
        raise Exception('url cannot be None')

    req = urllib2.Request(url)
    resp = urllib2.urlopen(req)

    last_modified = resp.info()['Last-Modified']
    logger.debug(url + ' last modified: ' + last_modified)
    lmdate = parse(last_modified)
    return is_datetime_after_celpp_week_start(lmdate)


def get_celpp_week_of_year_from_date(the_date):
    """ Returns CELPP week of year
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
        logger.debug(dir_to_create +
                     ' directory already exists skipping creation')
        return

    logger.debug('Creating ' + dir_to_create)
    os.makedirs(dir_to_create, mode=0o775)


def download_url_to_file(url, download_path, num_retries,
                         retry_sleep_time_secs):
    """Downloads `url` to path specified by `download_path`

    Downloads `url` with built in retry `num_retries` and sleep
    `retry_sleep_time_secs` between retries
    :param url: url to download file from
    :param download_path: full path to file to write output to
    :param num_retries: number of retries, 0 means none, 1 means one retry
    :param retry_sleep_time_secs: # seconds delay between retry, cannot
                                  be negative
    :raises: DownloadError if there was an error downloading the file or
             one of the parameters has invalid value
    """

    if url is None:
        raise DownloadError('url is not set')

    if download_path is None:
        raise DownloadError('download_path is not set')

    if num_retries is None:
        num_retries = 0

    if retry_sleep_time_secs is None:
        retry_sleep_time_secs = 0

    if retry_sleep_time_secs < 0:
        raise DownloadError('retry sleep time cannot be negative')

    count = 0
    while count <= num_retries:
        logger.debug('Try # ' + str(count) + ' of ' +
                     str(num_retries) + ' to download ' +
                     download_path + ' from ' + url)
        try:
            (f, info) = urllib.urlretrieve(url, download_path)
            return
        except Exception:
            logger.exception('Caught exception trying to download file')
        count += 1
        logger.debug('Download failed, sleeping ' +
                     str(retry_sleep_time_secs) +
                     'seconds')
        time.sleep(retry_sleep_time_secs)

    raise DownloadError('Unable to download file from ' +
                        url + ' to ' + download_path)


def gunzip_file(gzip_file, dest_file):
    """Uses Python gzip library to uncompress gzip file

       :param gzip_file: path to file to gunzip
       :param dest_file: gunzip destination file
       :raises IOError: If there was a problem with operation
    """
    f_in = gzip.open(gzip_file, 'rb')
    f_out = open(dest_file, 'wb')
    for f in f_in:
        f_out.write(f)
    f_out.flush()
    f_out.close()
    f_in.close()


def append_string_to_file(dest_file, the_string):
    """Appends text in `the_string` to file set by `dest_file`

       Opens file specified by `dest_file` in append mode and
       writes `the_string` contents to the file.  Function then
       flushes and closes file.

       :param dest_file: File to append data to
       :param the_string: Data to append to file
       :raises IOError: If there is an issue writing data to the file
       :raises TypeError: If `dest_file` is None
    """
    f = open(dest_file, 'a')
    f.write(the_string)
    f.flush()
    f.close()


def get_file_line_count(the_file):
    """Counts number of lines in file

       Opens file and counts number of lines in
       file
       :param the_file: file to examine
       :raises TypeError: if `the_file` is None
       :raises IOError: if there is an issue reading the file
    """
    f = open(the_file, 'r')
    counter = 0
    for line in f:
        counter += 1
    f.close()
    return counter


def run_external_command(cmd_to_run):
    """Runs command via external process
       Executes command in `cmd_to_run` which should
       be a single string command ie: ls -la or /bin/true
       :param cmd_to_run: The command with arguments to run
       :raises: All exceptions from subprocess.Popen()
       :returns: tuple (exitcode, stdout, stderr)
    """

    if cmd_to_run is None:
        return 256, '', 'Command must be set'

    logger.info("Running command " + cmd_to_run)

    p = subprocess.Popen(shlex.split(cmd_to_run),
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    out, err = p.communicate()
    return p.returncode, out, err


def setup_logging(theargs):
    """Sets up the logging for application

    Loggers are setup for:
    d3r.blastnfilter
    d3r.celpprunner
    d3r.celpp.blastnfilter
    d3r.celpp.dataimport
    d3r.celpp.glide
    d3r.celpp.makeblastdb
    d3r.celpp.proteinligprep
    d3r.celpp.evaluation
    d3r.celpp.task
    d3r.celpp.util
    d3r.celpp.uploader
    d3r.celpp.chimeraprep
    d3r.celpp.participant
    d3r.blast.ligand
    d3r.blast.hit
    d3r.blast.hit_sequence
    d3r.blast.query
    d3r.filter.filter
    d3r.utilities.analysis
    d3r.utilities.in_put
    d3r.utilities.run

    NOTE:  If new modules are added please add their loggers to this
    function

    The loglevel is set by theargs.loglevel and the format is set
    by LOG_FORMAT set at the top of this module.
    :param: theargs should have .loglevel set to one of the
    following strings: DEBUG, INFO, WARNING, ERROR, CRITICAL
    """
    theargs.logformat = LOG_FORMAT
    theargs.numericloglevel = logging.NOTSET
    if theargs.loglevel == 'DEBUG':
        theargs.numericloglevel = logging.DEBUG
    if theargs.loglevel == 'INFO':
        theargs.numericloglevel = logging.INFO
    if theargs.loglevel == 'WARNING':
        theargs.numericloglevel = logging.WARNING
    if theargs.loglevel == 'ERROR':
        theargs.numericloglevel = logging.ERROR
    if theargs.loglevel == 'CRITICAL':
        theargs.numericloglevel = logging.CRITICAL

    logger.setLevel(theargs.numericloglevel)
    logging.basicConfig(format=theargs.logformat)

    # There should be a line below for every package aka .py file
    # under celpp module
    logging.getLogger('d3r.blastnfilter')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpprunner')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celppreports')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.blastnfilter')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.challengedata')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.dataimport').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.glide').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.vina').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.makeblastdb')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.proteinligprep')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.evaluation').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.task').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.util').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.filetransfer')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.chimeraprep')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.celpp.participant')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.blast.ligand').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.blast.hit').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.blast.hit_sequence')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.blast.query').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.filter.filter').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.utilities.analysis')\
        .setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.utilities.in_put').setLevel(theargs.numericloglevel)
    logging.getLogger('d3r.utilities.run').setLevel(theargs.numericloglevel)
