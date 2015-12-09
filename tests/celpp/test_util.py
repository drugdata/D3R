#!/usr/bin/env python
# -*- coding: utf-8 -*-


import unittest
import tempfile
import os.path

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
            os.makedirs(os.path.join(temp_dir, '2015', 'dataset.week.1'), 0775)
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

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
