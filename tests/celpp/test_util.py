#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import os.path
import gzip

"""
test_task
----------------------------------

Tests for `task` module.
"""

import shutil
from datetime import date
from d3r.celpp import util
from d3r.celpp.util import DownloadError


class TestUtil(unittest.TestCase):
    def setUp(self):
        pass

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

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
