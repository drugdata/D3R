#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import os.path
import gzip
import stat
import logging

from urllib2 import URLError
from datetime import datetime
from dateutil.tz import tzutc
from dateutil.tz import tzlocal
from datetime import timedelta


"""
test_task
----------------------------------

Tests for `task` module.
"""

import shutil
from datetime import date
from d3r.celpp import util
from d3r.celpp.util import DownloadError
from d3r.celpp.task import D3RParameters


class TestUtil(unittest.TestCase):
    def setUp(self):
        pass

    def get_total_seconds(self, td):
        return (td.microseconds +
                (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

    def test_find_latest_year(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                fakefile = os.path.join(temp_dir, 'fakefile')
                open(fakefile, 'a').close()
                util.find_latest_year(fakefile)
                self.fail('Expected exception')
            except Exception:
                pass
            self.assertEqual(util.find_latest_year(temp_dir), None)

            os.mkdir(os.path.join(temp_dir, "foo"))
            os.mkdir(os.path.join(temp_dir, "2014"))
            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2014'))

            open(os.path.join(temp_dir, '2016'), 'a').close()
            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2014'))

            os.mkdir(os.path.join(temp_dir, '2012'))
            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2014'))

            os.mkdir(os.path.join(temp_dir, '2015'))

            self.assertEqual(util.find_latest_year(temp_dir),
                             os.path.join(temp_dir, '2015'))

        finally:
            shutil.rmtree(temp_dir)

    def test_find_latest_weekly_dataset(self):
        temp_dir = tempfile.mkdtemp()

        try:
            self.assertEqual(util.find_latest_weekly_dataset(temp_dir), None)

            os.mkdir(os.path.join(temp_dir, '2015'))

            self.assertEqual(util.find_latest_weekly_dataset(temp_dir), None)

            open(os.path.join(temp_dir, '2015', 'dataset.week.4'), 'a').close()

            self.assertEqual(util.find_latest_weekly_dataset(temp_dir), None)

            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.3'))

            self.assertEqual(util.find_latest_weekly_dataset(temp_dir),
                             os.path.join(temp_dir, '2015',
                                          'dataset.week.3'))

            # added this check to check fix for issue #13
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.9'))
            self.assertEqual(util.find_latest_weekly_dataset(temp_dir),
                             os.path.join(temp_dir, '2015',
                                          'dataset.week.9'))

            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week'))
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.'))

            # added this check to check fix for issue #13
            os.mkdir(os.path.join(temp_dir, '2015', 'dataset.week.10'))
            self.assertEqual(util.find_latest_weekly_dataset(temp_dir),
                             os.path.join(temp_dir, '2015',
                                          'dataset.week.10'))

        finally:
            shutil.rmtree(temp_dir)

    def test_get_celpp_week_of_year(self):
        # try passing none
        try:
            util.get_celpp_week_of_year_from_date(None)
            self.fail('Expected exception')
        except Exception:
            pass

        # try 1-1-2016
        celp_week = util.get_celpp_week_of_year_from_date(date(2016, 1, 1))
        self.assertEqual(celp_week[0], 1)
        self.assertEqual(celp_week[1], 2016)

        # try 12-18-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 12, 18))
        self.assertEqual(celp_week[0], 52)
        self.assertEqual(celp_week[1], 2015)

        # try 12-25-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 12, 25))
        self.assertEqual(celp_week[0], 53)
        self.assertEqual(celp_week[1], 2015)

        # try 12-31-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 12, 31))
        self.assertEqual(celp_week[0], 53)
        self.assertEqual(celp_week[1], 2015)

        # try 10-5-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 10, 5))
        self.assertEqual(celp_week[0], 41)
        self.assertEqual(celp_week[1], 2015)

        # try 10-8-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 10, 8))
        self.assertEqual(celp_week[0], 41)
        self.assertEqual(celp_week[1], 2015)

        # try 10-9-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 10, 9))
        self.assertEqual(celp_week[0], 42)
        self.assertEqual(celp_week[1], 2015)

        # try 10-10-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 10, 10))
        self.assertEqual(celp_week[0], 42)
        self.assertEqual(celp_week[1], 2015)

        # try 10-11-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 10, 11))
        self.assertEqual(celp_week[0], 42)
        self.assertEqual(celp_week[1], 2015)

        # try 10-12-2015
        celp_week = util.get_celpp_week_of_year_from_date(date(2015, 10, 12))
        self.assertEqual(celp_week[0], 42)
        self.assertEqual(celp_week[1], 2015)

    def test_create_celpp_week_dir(self):
        temp_dir = tempfile.mkdtemp()
        try:

            try:
                # try both parameters None
                util.create_celpp_week_dir(None, None)
                self.fail('Expected exception')
            except Exception:
                pass

            try:
                # try week tuple None
                util.create_celpp_week_dir(None, temp_dir)
            except Exception:
                pass

            try:
                # try celppdir None
                util.create_celpp_week_dir((1, 2015), None)
            except Exception:
                pass

            # try path already exists as dir
            os.makedirs(os.path.join(temp_dir, '2015', 'dataset.week.1'))
            util.create_celpp_week_dir((1, 2015), temp_dir)

            # try path already exists as file so it will fail
            open(os.path.join(temp_dir, '2015', 'dataset.week.2'), 'a').close()
            try:
                util.create_celpp_week_dir((2, 2015), temp_dir)
                self.fail('Expected OSError')
            except OSError:
                pass

            # try path successful create
            util.create_celpp_week_dir((3, 2015), temp_dir)
            new_dir = os.path.join(temp_dir, '2015', 'dataset.week.3')
            self.assertEquals(os.path.isdir(new_dir), True)
        finally:
            shutil.rmtree(temp_dir)

    def test_download_url_to_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # url is None
            try:
                util.download_url_to_file(None, None, 0, 0)
                self.fail('Expected DownloadError')
            except DownloadError as d:
                self.assertEquals(str(d), 'url is not set')

            # download_path is None
            try:
                util.download_url_to_file('foo', None, 0, 0)
                self.fail('Expected DownloadError')
            except DownloadError as d:
                self.assertEquals(str(d), 'download_path is not set')

            # retry sleep time negative
            try:
                util.download_url_to_file('foo', 'ha', 0, -1)
                self.fail('Expected DownloadError')
            except DownloadError as d:
                self.assertEquals(str(d), 'retry sleep time cannot be ' +
                                  'negative')

            # error download
            try:
                util.download_url_to_file('file://' + temp_dir +
                                          '/haha',
                                          os.path.join(temp_dir, 'boo'), 0, 0)
                self.fail('Expected DownloadError')
            except DownloadError as d:
                self.assertEquals(str(d), 'Unable to download file from ' +
                                  'file://' + temp_dir + '/haha to ' +
                                  os.path.join(temp_dir, 'boo'))

            # successful download retry sleep and num retries both None
            fake_file = os.path.join(temp_dir, 'hello')
            f = open(fake_file, 'w')
            f.write('hi\n')
            f.flush()
            f.close()
            util.download_url_to_file('file://' + fake_file,
                                      os.path.join(temp_dir, 'boo'),
                                      None, None)
            self.assertEquals(os.path.isfile(os.path.join(temp_dir, 'boo')),
                              True)

            # successful download retry sleep and num retries both 0
            util.download_url_to_file('file://' + fake_file,
                                      os.path.join(temp_dir, 'boo2'), 0, 0)
            self.assertEquals(os.path.isfile(os.path.join(temp_dir, 'boo2')),
                              True)

        finally:
            shutil.rmtree(temp_dir)

    def test_gunzip_file(self):
        temp_dir = tempfile.mkdtemp()
        try:

            # try invalid input
            try:
                util.gunzip_file(os.path.join(temp_dir, 'nonexistantfile'),
                                 os.path.join(temp_dir, 'ha'))
                self.fail('Expected IOError')
            except IOError:
                pass

            # try invalid output
            foogz = os.path.join(temp_dir, 'foo.gz')
            f = gzip.open(foogz, 'wb')
            f.write('hello\n')
            f.flush()
            f.close()
            try:

                util.gunzip_file(foogz,
                                 os.path.join(temp_dir, 'ha', 'ha'))
                self.fail('Expected IOError')
            except IOError:
                pass

            # try valid zip file
            util.gunzip_file(foogz,
                             os.path.join(temp_dir, 'foo'))

            f = open(os.path.join(temp_dir, 'foo'), 'r')
            self.assertEqual(f.readline(), 'hello\n')
            f.close()

        finally:
            shutil.rmtree(temp_dir)

    def test_append_string_to_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                # Try with double None
                util.append_string_to_file(None, None)
                self.fail('Expected TypeError')
            except TypeError:
                pass

            try:
                # Try with value set to None
                util.append_string_to_file(os.path.join(temp_dir, 'foo'),
                                           None)
            except TypeError:
                pass

            # try on non existant file
            util.append_string_to_file(os.path.join(temp_dir, 'non'), 'hi')
            f = open(os.path.join(temp_dir, 'non'))
            self.assertEqual(f.readline(), 'hi')
            f.close()

            # try on empty file
            open(os.path.join(temp_dir, 'empty'), 'a').close()
            util.append_string_to_file(os.path.join(temp_dir, 'empty'), 'hi')
            f = open(os.path.join(temp_dir, 'empty'))
            self.assertEqual(f.readline(), 'hi')
            f.close()

            # try on file with data
            util.append_string_to_file(os.path.join(temp_dir, 'empty'), 'how')
            f = open(os.path.join(temp_dir, 'empty'))
            self.assertEqual(f.readline(), 'hihow')

        finally:
            shutil.rmtree(temp_dir)

    def test_get_file_line_count(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                # Try with double None
                util.get_file_line_count(None)
                self.fail('Expected TypeError')
            except TypeError:
                pass

            try:
                # try on non existant file
                util.get_file_line_count(os.path.join(temp_dir, 'non'))
            except IOError:
                pass

            # try on empty file
            open(os.path.join(temp_dir, 'empty'), 'a').close()
            self.assertEqual(util.get_file_line_count(os.path.join(temp_dir,
                                                                   'empty')),
                             0)

            # try on file with data
            util.append_string_to_file(os.path.join(temp_dir, 'three'),
                                       'h\no\nw\n')
            self.assertEqual(util.get_file_line_count(os.path.join(temp_dir,
                                                                   'three')),
                             3)

        finally:
            shutil.rmtree(temp_dir)

    def test_get_all_celpp_years(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                util.get_all_celpp_years(os.path.join(temp_dir,
                                                      'doesnotexist'))
                self.fail('Expected exception')
            except Exception:
                pass

            # try with empty dir
            ylist = util.get_all_celpp_years(temp_dir)
            self.assertEqual(len(ylist), 0)

            # try on dir with no valid year directories
            os.mkdir(os.path.join(temp_dir, '123'))
            ylist = util.get_all_celpp_years(temp_dir)
            self.assertEqual(len(ylist), 0)

            # try on dir with 1 valid year
            os.mkdir(os.path.join(temp_dir, '2017'))
            ylist = util.get_all_celpp_years(temp_dir)
            self.assertEqual(len(ylist), 1)
            self.assertEqual(ylist[0], '2017')

            # try on dir with multiple valid year directories
            os.mkdir(os.path.join(temp_dir, '2016'))
            ylist = util.get_all_celpp_years(temp_dir)
            self.assertEqual(len(ylist), 2)
            ylist.sort()
            self.assertEqual(ylist[0], '2016')
            self.assertEqual(ylist[1], '2017')

        finally:
            shutil.rmtree(temp_dir)

    def test_get_celpp_year_from_path(self):

        # pass in None
        try:
            util.get_celpp_year_from_path(None)
            self.fail('Expected exception')
        except Exception:
            pass

        # pass in empty string
        self.assertEqual(util.get_celpp_year_from_path(''), '0')

        # pass in white space string
        self.assertEqual(util.get_celpp_year_from_path('  '), '0')

        # pass in path ending with year
        self.assertEqual(util.get_celpp_year_from_path('/home/foo/2013'),
                         '2013')

        # pass in path with invalid subdirectory on it
        self.assertEqual(util.get_celpp_year_from_path('/h/f/2013/cheese'),
                         '0')

        # pass in path with valid dataset week on end
        valid_week = os.path.join('h', '2013', util.DATA_SET_WEEK_PREFIX
                                  + '10')
        self.assertEqual(util.get_celpp_year_from_path(valid_week), '2013')

    def test_get_all_celpp_weeks(self):
        temp_dir = tempfile.mkdtemp()
        try:
            # try on invalid dir
            try:
                util.get_all_celpp_weeks(os.path.join(temp_dir,
                                                      'doesnotexist'))
                self.fail('Expected exception')
            except Exception:
                pass

            # try on empty dir
            week_list = util.get_all_celpp_weeks(temp_dir)
            self.assertEqual(len(week_list), 0)

            # try on dir with no valid weeks
            os.mkdir(os.path.join(temp_dir, 'foo'))
            week_list = util.get_all_celpp_weeks(temp_dir)
            self.assertEqual(len(week_list), 0)

            # try on dir with 1 valid week
            validweek = os.path.join(temp_dir, util.DATA_SET_WEEK_PREFIX
                                     + '12')
            os.mkdir(validweek)
            week_list = util.get_all_celpp_weeks(temp_dir)
            self.assertEqual(len(week_list), 1)
            self.assertEqual(week_list[0], '12')

            # try on dir with 2 valid weeks
            validweek = os.path.join(temp_dir, util.DATA_SET_WEEK_PREFIX
                                     + '51')
            os.mkdir(validweek)
            week_list = util.get_all_celpp_weeks(temp_dir)
            self.assertEqual(len(week_list), 2)
            week_list.sort()
            self.assertEqual(week_list[0], '12')
            self.assertEqual(week_list[1], '51')

        finally:
            shutil.rmtree(temp_dir)

    def test_get_celpp_week_number_from_path(self):
        try:
            util.get_celpp_week_number_from_path(None)
            self.fail('Expected exception')
        except Exception:
            pass

        # try with empty string
        foo = util.get_celpp_week_number_from_path('')
        self.assertEqual(foo, '0')

        # try with path missing dataset.week. prefix
        foo = util.get_celpp_week_number_from_path('/h/b/hi')
        self.assertEqual(foo, '0')

        # try with valid path
        foo = util.get_celpp_week_number_from_path('/h/b/'
                                                   + util.DATA_SET_WEEK_PREFIX
                                                   + '7')
        self.assertEqual(foo, '7')

        foo = util.get_celpp_week_number_from_path('/h/b/'
                                                   + util.DATA_SET_WEEK_PREFIX
                                                   + '12')
        self.assertEqual(foo, '12')

    def test_get_previous_date_of_previous_friday_from_date(self):
        # test None passed in
        try:
            util.get_previous_friday_from_date(None)
            self.fail('Expected exception')
        except Exception:
            pass

        # test non date passed in
        try:
            util.get_previous_friday_from_date('hello')
            self.fail('Expected exception')
        except Exception:
            pass

        # test date on Friday
        dt = datetime(2016, 6, 17, 8, 0, 5)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2016)
        self.assertTrue(newdt.day == 17)
        self.assertTrue(newdt.month == 6)
        self.assertTrue(newdt.hour == 8)
        self.assertTrue(newdt.minute == 0)
        self.assertTrue(newdt.second == 5)

        # test date on Saturday
        dt = datetime(2015, 8, 1, 15, 20, 2)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2015)
        self.assertTrue(newdt.day == 31)
        self.assertTrue(newdt.month == 7)
        self.assertTrue(newdt.hour == 15)
        self.assertTrue(newdt.minute == 20)
        self.assertTrue(newdt.second == 2)

        # test date on Sunday
        dt = datetime(2016, 6, 19, 8, 0, 5)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2016)
        self.assertTrue(newdt.day == 17)
        self.assertTrue(newdt.month == 6)

        # test date on Monday
        dt = datetime(2016, 6, 20, 8, 0, 5)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2016)
        self.assertTrue(newdt.day == 17)
        self.assertTrue(newdt.month == 6)

        # test date on Tuesday
        dt = datetime(2016, 6, 21, 8, 0, 5)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2016)
        self.assertTrue(newdt.day == 17)
        self.assertTrue(newdt.month == 6)

        # test date on Wednesday
        dt = datetime(2016, 6, 22, 8, 0, 5)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2016)
        self.assertTrue(newdt.day == 17)
        self.assertTrue(newdt.month == 6)

        # test date on Thursday
        dt = datetime(2016, 6, 23, 8, 0, 5)
        newdt = util.get_previous_friday_from_date(dt)
        self.assertTrue(newdt.year == 2016)
        self.assertTrue(newdt.day == 17)
        self.assertTrue(newdt.month == 6)

    def test_is_datetime_after_celpp_week_start(self):

        try:
            util.is_datetime_after_celpp_week_start(None)
            self.fail('Expected exception')
        except Exception:
            pass
        try:
            util.is_datetime_after_celpp_week_start('hi')
            self.fail('Expected exception')
        except AssertionError:
            pass

        prev_friday = util.get_previous_friday_from_date(
            datetime.now(tzlocal()))
        oneday = timedelta(days=1)
        sat = prev_friday + oneday
        self.assertTrue(util.is_datetime_after_celpp_week_start(sat))
        thurs = prev_friday - oneday
        self.assertFalse(util.is_datetime_after_celpp_week_start(thurs))

    def test_has_url_been_updated_since_start_of_celpp_week(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                util.has_url_been_updated_since_start_of_celpp_week(None)
                self.fail('Expected exception')
            except Exception:
                pass

            fakefile = os.path.join(temp_dir, 'foo')

            try:
                util.has_url_been_updated_since_start_of_celpp_week('file://' +
                                                                    fakefile)
                self.fail('Expected exception')
            except URLError:
                pass

            f = open(fakefile, 'w')
            f.write('hi\n')
            f.flush()
            f.close()

            prev_friday = util.get_previous_friday_from_date(
                datetime.now(tzlocal()))

            thurs = prev_friday - timedelta(days=1)
            dse = thurs - datetime(1970, 01, 01, tzinfo=tzutc())
            secs_since_epoch = self.get_total_seconds(dse)
            os.utime(fakefile, (secs_since_epoch, secs_since_epoch))
            val = util.has_url_been_updated_since_start_of_celpp_week(
                'file://' + fakefile)
            self.assertEqual(val, False)

            sat = prev_friday + timedelta(days=1)
            dse = sat - datetime(1970, 01, 01, tzinfo=tzutc())
            secs_since_epoch = self.get_total_seconds(dse)
            os.utime(fakefile, (secs_since_epoch, secs_since_epoch))
            val = util.has_url_been_updated_since_start_of_celpp_week(
                'file://' + fakefile)
            self.assertEqual(val, True)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_command_not_set(self):
        ecode, out, err = util.run_external_command(None)
        self.assertEqual(ecode, 256)
        self.assertEqual(out, '')
        self.assertEqual(err, 'Command must be set')

    def test_run_external_command_cmd_does_not_exist(self):
        temp_dir = tempfile.mkdtemp()
        try:
            try:
                noexist = os.path.join(temp_dir, 'noexist')
                ecode, out, err = util.run_external_command(noexist)
                self.fail('Expected OSError')
            except OSError:
                pass

        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_success_with_output(self):
        temp_dir = tempfile.mkdtemp()
        try:
            script = os.path.join(temp_dir, 'yo.py')

            # create a small python script that outputs args passed
            # in to standard out, writes error to standard error
            #  and exits with 0 exit code
            f = open(script, 'w')
            f.write('#! /usr/bin/env python\n\n')
            f.write('import sys\n')
            f.write('sys.stdout.write(sys.argv[1])\n')
            f.write('sys.stdout.write(sys.argv[2])\n')
            f.write('sys.stderr.write("error")\n')
            f.write('sys.exit(0)\n')
            f.flush()
            f.close()
            os.chmod(script, stat.S_IRWXU)

            ecode, out, err = util.run_external_command(script + ' hi how')

            self.assertEqual(err, 'error')
            self.assertEqual(out, 'hihow')
            self.assertEqual(ecode, 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_fail_with_output(self):
        temp_dir = tempfile.mkdtemp()
        try:
            script = os.path.join(temp_dir, 'yo.py')

            # create a small python script that outputs args passed
            # in to standard out, writes error to standard error
            #  and exits with 0 exit code
            f = open(script, 'w')
            f.write('#! /usr/bin/env python\n\n')
            f.write('import sys\n')
            f.write('sys.stdout.write(sys.argv[1])\n')
            f.write('sys.stdout.write(sys.argv[2])\n')
            f.write('sys.stderr.write("2error")\n')
            f.write('sys.exit(2)\n')
            f.flush()
            f.close()
            os.chmod(script, stat.S_IRWXU)

            ecode, out, err = util.run_external_command(script + ' hi how')

            self.assertEqual(err, '2error')
            self.assertEqual(out, 'hihow')
            self.assertEqual(ecode, 2)

        finally:
            shutil.rmtree(temp_dir)

    def test_setup_logging(self):
        logger = logging.getLogger('funlogger')
        theargs = D3RParameters()
        theargs.loglevel = 'INFO'
        util.setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.task')
                         .getEffectiveLevel(),
                         logging.INFO)
        self.assertEqual(theargs.numericloglevel, logging.INFO)
        logger.debug('test')

        theargs.loglevel = 'DEBUG'
        util.setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.task')
                         .getEffectiveLevel(),
                         logging.DEBUG)
        self.assertEqual(theargs.numericloglevel, logging.DEBUG)

        theargs.loglevel = 'WARNING'
        util.setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.task')
                         .getEffectiveLevel(),
                         logging.WARNING)
        self.assertEqual(theargs.numericloglevel, logging.WARNING)

        theargs.loglevel = 'ERROR'
        util.setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.task')
                         .getEffectiveLevel(),
                         logging.ERROR)
        self.assertEqual(theargs.numericloglevel, logging.ERROR)

        theargs.loglevel = 'CRITICAL'
        util.setup_logging(theargs)
        self.assertEqual(logging.getLogger('d3r.celpp.task')
                         .getEffectiveLevel(),
                         logging.CRITICAL)
        self.assertEqual(theargs.numericloglevel, logging.CRITICAL)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
