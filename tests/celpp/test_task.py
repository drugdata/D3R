#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_task
----------------------------------

Tests for `task` module.
"""
import unittest
import tempfile
import shutil
import platform
import os
import pwd
import gzip
from email.mime.multipart import MIMEMultipart
import configparser
from mock import Mock

from d3r.celpp.task import D3RParameters
from d3r.celpp.task import UnsetPathError
from d3r.celpp.task import UnsetStageError
from d3r.celpp.task import UnsetNameError
from d3r.celpp.task import UnsetCommandError
from d3r.celpp.task import UnsetFileNameError
from d3r.celpp.task import D3RTask
from d3r.celpp.task import SmtpEmailer
from d3r.celpp.task import EmailSendError
from d3r.celpp.task import Attachment
from d3r.celpp.filetransfer import FtpFileTransfer
from d3r.celpp.task import SmtpConfig
from d3r.celpp.task import SmtpEmailerFactory
from d3r.celpp.task import WebsiteServiceConfig


class MockException(Exception):
    pass


class MockFileTransfer(FtpFileTransfer):
    """Fake FileUploader
    """
    def __init__(self):
        self._upload_files_except_object = None
        self._upload_summary_except_object = None
        self._upload_summary = None
        self._upload_files = None
        self._connect = True

    def set_upload_files(self, val):
        self._upload_files = val

    def set_upload_files_exception(self, except_object):
        self._upload_files_except_object = except_object

    def upload_files(self, list_of_files):
        """

        :param list_of_files:
        :return:
        """
        if self._upload_files_except_object is not None:
            raise self._upload_files_except_object
        return self._upload_files

    def set_upload_summary(self, val):
        self._upload_summary = val

    def set_upload_summary_exception(self, except_object):
        self._upload_summary_except_object = except_object

    def get_upload_summary(self):
        """

        :return: String set in set_upload_summary()
        """
        if self._upload_summary_except_object is not None:
            raise self._upload_summary_except_object

        return self._upload_summary

    def set_connect(self, val):
        self._connect = val

    def connect(self):
        return self._connect

    def disconnect(self):
        pass


class TestD3rTask(unittest.TestCase):

    def setUp(self):
        pass

    def test_constructor(self):
        params = D3RParameters()
        task = D3RTask('/path', params)
        self.assertEqual(task.get_name(), None)
        self.assertEqual(task.get_path(), '/path')
        self.assertEqual(task.get_stage(), None)
        self.assertEqual(task.get_status(), D3RTask.UNKNOWN_STATUS)
        self.assertEqual(task.get_error(), None)
        task.set_name('foo')
        task.set_path('blah')
        task.set_stage(4)
        task.set_status(D3RTask.START_STATUS)
        task.set_error('error')
        task.set_file_transfer('yah')

        self.assertEqual(task.get_name(), 'foo')
        self.assertEqual(task.get_path(), 'blah')
        self.assertEqual(task.get_stage(), 4)
        self.assertEqual(task.get_status(), D3RTask.START_STATUS)
        self.assertEqual(task.get_error(), 'error')
        self.assertEqual(task.get_file_transfer(), 'yah')

        params.ftpconfig = '/somefile'
        task = D3RTask('/path', params)
        self.assertEqual(task.get_file_transfer(), None)

    def test_get_dir_name(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.get_dir_name()
            self.fail('Expected UnsetStageError')
        except UnsetStageError:
            pass

        task.set_stage(1)
        try:
            task.get_dir_name()
            self.fail('Expected UnsetNameError')
        except UnsetNameError:
            pass

        task.set_name('foo')
        self.assertEqual(task.get_dir_name(), 'stage.1.foo')

    def test_get_dir(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.get_dir()
            self.fail('Expected UnsetPathError')
        except UnsetPathError:
            pass

        task.set_path('/blah')

        try:
            task.get_dir()
            self.fail('Expected UnsetStageError')
        except UnsetStageError:
            pass

        task.set_stage(1)
        try:
            task.get_dir()
            self.fail('Expected UnsetNameError')
        except UnsetNameError:
            pass

        task.set_name('foo')

        self.assertEqual(task.get_dir(), '/blah/stage.1.foo')

    def test_get_set_email_log(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        self.assertEqual(task.get_email_log(), None)
        task.set_email_log('hi')
        self.assertEqual(task.get_email_log(), 'hi')

    def test_get_program_name(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        self.assertTrue(task._get_program_name().find('task.py') > 0,
                        task._get_program_name())
        params.program = 'proggy'
        params.version = 'versy'
        task = D3RTask(None, params)
        self.assertEqual(task._get_program_name(), 'proggy versy')

    def test_get_program_version(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        self.assertEqual(task.get_program_version(), '')
        params.version = '0.1.1'
        task = D3RTask(None, params)
        self.assertEqual(task.get_program_version(), '0.1.1')

    def test_get_uploadable_files(self):
        task = D3RTask(None, D3RParameters())
        self.assertEqual(task.get_uploadable_files(), [])
        temp_dir = tempfile.mkdtemp()
        try:
            task = D3RTask(temp_dir, D3RParameters())
            task.set_stage(1)
            task.set_name('foo')
            task.create_dir()
            self.assertEqual(task.get_uploadable_files(), [])

            # add error file
            err_file = os.path.join(task.get_dir(), D3RTask.ERROR_FILE)
            open(err_file, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 1)
            self.assertEqual(flist[0], err_file)

            # add stderr file
            foostderr = os.path.join(task.get_dir(), 'foo.stderr')
            open(foostderr, 'a').close()
            flist = task.get_uploadable_files()
            flist.index(err_file)
            flist.index(foostderr)
            self.assertEqual(len(flist), 2)

            # add stdout file
            foostdout = os.path.join(task.get_dir(), 'foo.stdout')
            open(foostdout, 'a').close()
            flist = task.get_uploadable_files()
            self.assertEqual(len(flist), 3)
            flist.index(err_file)
            flist.index(foostderr)
            flist.index(foostdout)
        finally:
            shutil.rmtree(temp_dir)

    def test_can_run(self):
        task = D3RTask(None, D3RParameters())
        self.assertEqual(task.can_run(), False)

    def test_upload_task(self):
        task = D3RTask(None, D3RParameters())
        self.assertEqual(task._file_uploader, None)
        # try calling upload_task where _file_uploader is None
        task._upload_task()

        # try calling upload_task where _file_uploader is NOT None
        params = D3RParameters()
        params.ftpconfig = 'foo'
        task = D3RTask(None, params)
        self.assertEqual(task._file_uploader, None)
        # try calling upload_task where _file_uploader is None
        task._upload_task()

        # try calling upload_task where _file_uploader.upload_files
        # throws exception
        uploader = MockFileTransfer()
        uploader.set_upload_files_exception(MockException('hi'))
        task = D3RTask(None, D3RParameters())
        task.set_file_transfer(uploader)
        task._upload_task()

        # try calling where upload_summary throws exception
        uploader.set_upload_files_exception(None)
        uploader.set_upload_summary_exception(MockException('boo'))
        task._upload_task()

        # try calling where upload_summary returns None
        uploader.set_upload_summary_exception(None)
        uploader.set_upload_summary(None)
        task._upload_task()

        uploader.set_upload_summary('hello')
        task._upload_task()
        self.assertEqual(task._email_log, '\nhello\n')

    def test_run(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        self.assertEqual(task._can_run, None)
        task.run()
        self.assertEqual(task._can_run, False)
        task._can_run = True

        temp_dir = tempfile.mkdtemp()
        task.set_name('foo')
        task.set_stage(1)
        task.set_path(temp_dir)
        try:
            task.run()
        finally:
            shutil.rmtree(temp_dir)

    def test_write_to_file(self):
        tempDir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            try:
                task.write_to_file('hello', None)
                self.fail('Expected UnsetFileNameError')
            except UnsetFileNameError:
                pass
            try:
                task.write_to_file('hello', 'foo')
                self.fail('Expected UnsetPathError')
            except UnsetPathError:
                pass
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(tempDir)

            try:
                task.write_to_file('hello', 'foo')
                self.fail('Expected IOError')
            except IOError:
                pass
            task.create_dir()
            task.write_to_file('hello', 'foo')
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                            'foo')), True)
        finally:
            shutil.rmtree(tempDir)

    def test_create_dir(self):
        tempDir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            try:
                task.create_dir()
                self.fail('Expected UnsetPathError')
            except UnsetPathError:
                pass
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(tempDir)
            self.assertEqual(task.create_dir(),
                             os.path.join(tempDir, 'stage.1.foo'))
        finally:
            shutil.rmtree(tempDir)

    def test_update_status_from_filesystem(self):
        params = D3RParameters()
        task = D3RTask(None, params)

        try:
            task.update_status_from_filesystem()
            self.fail("Expected UnsetPathError")
        except UnsetPathError:
            pass

        tempDir = tempfile.mkdtemp()
        try:
            task.set_path(tempDir)

            # Test unset stage
            try:
                task.update_status_from_filesystem()
                self.fail("Expected UnsetStageError")
            except UnsetStageError:
                pass
        finally:
            shutil.rmtree(tempDir)

        task.set_stage(1)
        task.set_name('foo')

        try_update_status_from_filesystem(self, task)

    def test_start_no_version(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(temp_dir, params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                            D3RTask.START_FILE)), True)
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.START_STATUS)
            task.start()
            self.assertNotEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask
                                                         .ERROR_FILE)), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_start_with_version(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.version = '1.6.0'
            task = D3RTask(temp_dir, params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
            sfile = os.path.join(task.get_dir(),
                                 D3RTask.START_FILE)

            self.assertEqual(os.path.isfile(sfile), True)
            f = open(sfile, 'r')
            self.assertEqual(f.read(), '1.6.0')
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.START_STATUS)
            task.start()
            self.assertNotEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask
                                                         .ERROR_FILE)), True)

        finally:
            shutil.rmtree(temp_dir)

    def test_end(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(temp_dir, params)
            task.set_stage(1)
            task.set_name('foo')
            task.start()
            task.end()
            self.assertEqual(os.path.isfile(os.path.join(task.get_dir(),
                                                         D3RTask.
                                                         COMPLETE_FILE)), True)
            self.assertEqual(task.get_error(), None)
            self.assertEqual(task.get_status(), D3RTask.COMPLETE_STATUS)
            task.set_error('some error')
            task.end()
            self.assertEqual(task.get_error(), 'some error')
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)

            err_file = os.path.join(task.get_dir(), D3RTask.ERROR_FILE)

            self.assertEqual(os.path.isfile(err_file), True)
            f = open(err_file, 'r')
            err_msg = f.readline()
            self.assertEqual(err_msg, 'some error\n')
            task.set_status(D3RTask.ERROR_STATUS)
            os.remove(err_file)
            task.set_error(None)
            task.end()
            self.assertEqual(task.get_status(), D3RTask.ERROR_STATUS)
            self.assertEqual(os.path.isfile(err_file), True)
            self.assertEqual(os.path.getsize(err_file), 0)

        finally:
            shutil.rmtree(temp_dir)

    def test_send_email_where_email_is_none(self):
        params = D3RParameters()
        params.email = None
        task = D3RTask('/foo', params)
        task.set_stage(1)
        task.set_name('foo')
        task._send_end_email()

    def test_send_end_email_throws_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            params.smtp = 'doesnotexistxxx.x.x'
            params.smtpport = 12345
            params.email = 'foo@foo.com'
            task = D3RTask(temp_dir, params)
            task.set_stage(1)
            task.set_name('foo')
            task._send_end_email()

        finally:
            shutil.rmtree(temp_dir)

    def test_get_mail_truncated_string_val_is_none(self):
        params = D3RParameters()
        task = D3RTask('/foo', params)
        self.assertEqual(task._get_email_truncated_string(None,
                                                          100), None)

    def test_get_mail_truncated_string_maxchars_none_or_negative(self):
        params = D3RParameters()
        task = D3RTask('/foo', params)
        self.assertEqual(task._get_email_truncated_string('foo',
                                                          None), 'foo')

        self.assertEqual(task._get_email_truncated_string('foo',
                                                          -1), 'foo')

    def test_get_mail_truncated_string_truncated(self):
        params = D3RParameters()
        task = D3RTask('/foo', params)
        val = '123456789 abcdefghijklmnop'
        res = task._get_email_truncated_string(val, 15)
        self.assertEqual(res, D3RTask.TEXT_TRUNCATED_STR +
                         'bcdefghijklmnop')

    def test_run_external_command_all_params_None(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.run_external_command(None, None, None)
            self.fail('expected UnsetNameError')
        except UnsetNameError:
            pass

    def test_run_external_command_name_None(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.run_external_command(None, 'true', None)
            self.fail('expected UnsetNameError')
        except UnsetNameError:
            pass

    def test_run_external_command_name_cmd_to_run_None(self):
        params = D3RParameters()
        task = D3RTask(None, params)
        try:
            task.run_external_command('hello', None, None)
            self.fail('expected UnsetCommandError')
        except UnsetCommandError:
            pass

    def test_run_external_command_name_failure_is_fatal_is_None(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(temp_dir)
            task.create_dir()
            task.set_path(temp_dir)
            self.assertEquals(1, task.run_external_command('hi',
                                                           'false',
                                                           None))
            self.assertEquals(task.get_error(), 'Non zero exit code:' +
                              ' 1 received. Standard out:  Standard error: ')
            self.assertEquals(task.get_status(), D3RTask.ERROR_STATUS)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_cmd_raises_exception(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(temp_dir)
            task.create_dir()
            task.set_path(temp_dir)
            self.assertEquals(1, task.run_external_command('hi',
                                                           'asdfasdf',
                                                           None))
            self.assertEquals(task.get_error(), 'Caught Exception trying ' +
                              'to run asdfasdf : [Errno 2] No such file ' +
                              'or directory')

            self.assertEquals(os.path.exists(os.path.join(task.get_dir(),
                                                          'hi.stdout')),
                              False)
            self.assertEquals(os.path.exists(os.path.join(task.get_dir(),
                                                          'hi.stderr')),
                              False)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_cmd_fails_and_is_fatal(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(temp_dir)
            task.create_dir()
            task.set_path(temp_dir)
            self.assertEquals(1, task.run_external_command('hi',
                                                           'false',
                                                           True))
            self.assertEquals(task.get_error(), 'Non zero exit code:' +
                              ' 1 received. Standard out:  Standard error: ')
            self.assertEquals(task.get_status(), D3RTask.ERROR_STATUS)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_cmd_fails_and_isnot_fatal(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(temp_dir)
            task.create_dir()
            task.set_path(temp_dir)
            self.assertEquals(1, task.run_external_command('hi',
                                                           'false',
                                                           False))
            self.assertEquals(task.get_error(), None)
            self.assertEquals(task.get_status(), D3RTask.UNKNOWN_STATUS)
        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_cmd_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(temp_dir)
            task.create_dir()
            task.set_path(temp_dir)
            self.assertEquals(0, task.run_external_command('hi',
                                                           'echo hi',
                                                           False))
            self.assertEquals(task.get_error(), None)
            self.assertEquals(task.get_status(), D3RTask.UNKNOWN_STATUS)
            self.assertEquals(os.path.exists(os.path.join(task.get_dir(),
                                                          'hi.stdout')),
                              True)
            self.assertEquals(os.path.exists(os.path.join(task.get_dir(),
                                                          'hi.stderr')),
                              True)

        finally:
            shutil.rmtree(temp_dir)

    def test_run_external_command_withtimeout_cmd_succeeds(self):
        temp_dir = tempfile.mkdtemp()
        try:
            params = D3RParameters()
            task = D3RTask(None, params)
            task.set_name('foo')
            task.set_stage(1)
            task.set_path(temp_dir)
            task.create_dir()
            task.set_path(temp_dir)
            self.assertEquals(0, task.run_external_command('hi',
                                                           'echo hi',
                                                           False,
                                                           timeout=30))
            self.assertEquals(task.get_error(), None)
            self.assertEquals(task.get_status(), D3RTask.UNKNOWN_STATUS)
            self.assertEquals(os.path.exists(os.path.join(task.get_dir(),
                                                          'hi.stdout')),
                              True)
            self.assertEquals(os.path.exists(os.path.join(task.get_dir(),
                                                          'hi.stderr')),
                              True)

        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_emailer_generate_from_address_using_login_and_host(self):
        emailer = SmtpEmailer()
        val = emailer._generate_from_address_using_login_and_host()
        self.assertEqual(val, pwd.getpwuid(os.getuid())[0] + '@' +
                         platform.node())

        val = emailer._generate_from_address_using_login_and_host(hostname='')
        self.assertEqual(val, pwd.getpwuid(os.getuid())[0] + '@localhost')

    def test_sending_real_email_but_to_invalid_port_and_host(self):
        emailer = SmtpEmailer(port=1231231231232)
        try:
            emailer.send_email('bob@bob.com', ['joe@joe.com'], 'subby2',
                               'hi\n')
            self.fail('Expected EmailSendError')
        except EmailSendError:
            pass

    def test_get_server_with_altserver_set(self):
        emailer = SmtpEmailer()
        fake = D3RParameters()
        fake.hi = 'fake'
        emailer.set_alternate_smtp_server(fake)
        self.assertEqual(emailer._get_server().hi, fake.hi)

    def test_get_server_with_altserver_set_and_login_needed(self):
        emailer = SmtpEmailer(user='bob', password='smith')
        mockserver = D3RParameters()
        mockserver.login = Mock()
        emailer.set_alternate_smtp_server(mockserver)
        emailer._get_server()
        mockserver.login.assert_called_once_with('bob', 'smith')

    def test_send_valid_email_to_fake_server(self):
        emailer = SmtpEmailer()
        mockserver = D3RParameters()
        mockserver.sendmail = Mock()
        mockserver.quit = Mock()

        emailer.set_alternate_smtp_server(mockserver)
        emailer.send_email(['joe@joe.com'], 'subby2', 'hi\n')
        mockserver.quit.assert_any_call()
        self.assertEqual(mockserver.sendmail.call_count, 1)

    def test_send_invalid_email_to_fake_server_that_throws_exception(self):
        emailer = SmtpEmailer()
        mockserver = D3RParameters()
        mockserver.sendmail = Mock(side_effect=IOError('some error'))
        mockserver.quit = Mock()

        emailer.set_alternate_smtp_server(mockserver)
        try:
            emailer.send_email(['joe@joe.com'], 'subby2', 'hi\n')
        except EmailSendError:
            pass

        mockserver.quit.assert_any_call()

    def test_attachment_class(self):
        a = Attachment('/foo', 'hi')
        self.assertEqual(a.get_desired_name(), 'hi')
        self.assertEqual(a.get_file_to_attach(), '/foo')

    def test_append_attachments(self):
        temp_dir = tempfile.mkdtemp()
        try:
            emailer = SmtpEmailer()
            # try no attachments
            msg_root = MIMEMultipart("alternative")
            res = emailer._append_attachments(msg_root, None)
            self.assertEqual(msg_root, res)

            # try 1 text attachment
            txtfile = os.path.join(temp_dir, 'my.txt')
            f = open(txtfile, 'w')
            f.write('hello')
            f.flush()
            f.close()
            msg_root = MIMEMultipart("alternative")
            a = Attachment(txtfile, 'renamed.txt')
            res = emailer._append_attachments(msg_root, [a])
            res.as_string().index('Content-Type: text/plain; '
                                  'charset="us-ascii"')
            res.as_string().index('Content-Disposition: attachment; '
                                  'filename="renamed.txt"')
            res.as_string().index('hello')

            # try 1 image attachment
            img = os.path.join(temp_dir, 'yo.ppm')
            f = open(img, 'w')
            f.write('P6\n1 1\n255\n')
            f.write('%c' % 255)
            f.flush()
            f.close()
            msg_root = MIMEMultipart("alternative")
            b = Attachment(img, None)
            res = emailer._append_attachments(msg_root, [b])
            res.as_string().index('Content-Type: image/x-portable-pixmap')
            res.as_string().index('Content-Disposition: attachment; '
                                  'filename="yo.ppm"')

            # try 1 gzip attachment
            mygz = os.path.join(temp_dir, 'foo.gz')
            f = gzip.open(mygz, 'wb')
            f.write('some compressed data')
            f.flush()
            f.close()
            c = Attachment(mygz, 'well.gz')
            res = emailer._append_attachments(msg_root, [c])
            res.as_string().index('Content-Type: application/octet-stream')
            res.as_string().index('Content-Disposition: attachment; '
                                  'filename="well.gz"')

            # try 3 attachments above together
            msg_root = MIMEMultipart("alternative")
            res = emailer._append_attachments(msg_root, [a, b, c])
            res.as_string().index('Content-Type: text/plain; '
                                  'charset="us-ascii"')
            res.as_string().index('Content-Disposition: attachment; '
                                  'filename="renamed.txt"')
            res.as_string().index('hello')

            res.as_string().index('Content-Type: image/x-portable-pixmap')
            res.as_string().index('Content-Disposition: attachment; '
                                  'filename="yo.ppm"')

            res.as_string().index('Content-Type: application/octet-stream')
            res.as_string().index('Content-Disposition: attachment; '
                                  'filename="well.gz"')
        finally:
            shutil.rmtree(temp_dir)

    def test_websiteserviceconfig_noconfig(self):
        webconfig = WebsiteServiceConfig()
        self.assertEqual(webconfig.get_timeout(), 0.1)
        self.assertEqual(webconfig.get_rmsd_url(), None)
        self.assertEqual(webconfig.get_basicauth_password(), None)
        self.assertEqual(webconfig.get_basicauth_user(), None)
        self.assertEqual(webconfig.get_apikey(), None)
        self.assertEqual(webconfig.get_portal_name(), 'notset')
        self.assertEqual(webconfig.get_source(), 'notset')
        self.assertEqual(webconfig.get_targets_url(), None)

    def test_websiteserviceconfig_nonexistantconfig(self):
        temp_dir = tempfile.mkdtemp()
        try:
            nonexist = os.path.join(temp_dir, 'foo')
            webconfig = WebsiteServiceConfig(configfile=nonexist)
            self.assertEqual(webconfig.get_timeout(), 0.1)
            self.assertEqual(webconfig.get_rmsd_url(), None)
            self.assertEqual(webconfig.get_basicauth_password(), None)
            self.assertEqual(webconfig.get_basicauth_user(), None)
            self.assertEqual(webconfig.get_apikey(), None)
            self.assertEqual(webconfig.get_portal_name(), 'notset')
            self.assertEqual(webconfig.get_source(), 'notset')
            self.assertEqual(webconfig.get_targets_url(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_websiteserviceconfig_emptyconfig(self):
        temp_dir = tempfile.mkdtemp()
        try:
            empty = os.path.join(temp_dir, 'foo')
            open(empty, 'a').close()
            webconfig = WebsiteServiceConfig(configfile=empty)
            self.assertEqual(webconfig.get_timeout(), 0.1)
            self.assertEqual(webconfig.get_rmsd_url(), None)
            self.assertEqual(webconfig.get_basicauth_password(), None)
            self.assertEqual(webconfig.get_basicauth_user(), None)
            self.assertEqual(webconfig.get_apikey(), None)
            self.assertEqual(webconfig.get_portal_name(), None)
            self.assertEqual(webconfig.get_source(), None)
            self.assertEqual(webconfig.get_targets_url(), None)

        finally:
            shutil.rmtree(temp_dir)

    def test_websiteserviceconfig_validconfig(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(WebsiteServiceConfig.DEFAULT)
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_URL,
                    'http://localhost')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_TIMEOUT,
                    '25.1')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_PORTAL_NAME,
                    'portalname')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_SOURCE,
                    'source')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_BASIC_PASS,
                    'thepass')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_BASIC_USER,
                    'user')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_APIKEY,
                    'key')
            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()

            webconfig = WebsiteServiceConfig(configfile=cfile)
            self.assertEqual(webconfig.get_timeout(), 25.1)
            self.assertEqual(webconfig.get_rmsd_url(),
                             'http://localhost/rmsd')
            self.assertEqual(webconfig.get_basicauth_password(), 'thepass')
            self.assertEqual(webconfig.get_basicauth_user(), 'user')
            self.assertEqual(webconfig.get_apikey(), 'key')
            self.assertEqual(webconfig.get_portal_name(), 'portalname')
            self.assertEqual(webconfig.get_source(), 'source')
            self.assertEqual(webconfig.get_targets_url(),
                             'http://localhost/week')
        finally:
            shutil.rmtree(temp_dir)

    def test_websiteserviceconfig_urlwithslashatendandinvalidtimeout(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(WebsiteServiceConfig.DEFAULT)
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_URL,
                    'http://localhost/')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_TIMEOUT,
                    'hi')

            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()

            webconfig = WebsiteServiceConfig(configfile=cfile)
            self.assertEqual(webconfig.get_timeout(), 0.1)
            self.assertEqual(webconfig.get_rmsd_url(),
                             'http://localhost/rmsd')
            self.assertEqual(webconfig.get_targets_url(),
                             'http://localhost/week')
        finally:
            shutil.rmtree(temp_dir)

    def test_websiteserviceconfig_parse_error(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(WebsiteServiceConfig.DEFAULT)
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_URL,
                    'http://localhost/')
            con.set(WebsiteServiceConfig.DEFAULT,
                    WebsiteServiceConfig.WEB_TIMEOUT,
                    'hi')

            with open(cfile, 'w') as f:
                f.write('[smtp\nas = xx=12=\n')
                f.flush()
            webconfig = WebsiteServiceConfig(configfile=cfile)
            self.assertEqual(webconfig.get_timeout(), 0.1)
        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_config_noconfig(self):
        smtpconfig = SmtpConfig()
        self.assertEqual(smtpconfig.get_from_address(), None)
        self.assertEqual(smtpconfig.get_host(), 'localhost')
        self.assertEqual(smtpconfig.get_port(), 25)
        self.assertEqual(smtpconfig.get_password(), None)
        self.assertEqual(smtpconfig.get_from_address(), None)
        self.assertEqual(smtpconfig.get_replyto_address(), None)

    def test_smtp_config_nonexistantconfigfile(self):
        temp_dir = tempfile.mkdtemp()
        try:
            nonexist_file = os.path.join(temp_dir, 'foo')
            smtpconfig = SmtpConfig(configfile=nonexist_file)
            self.assertEqual(smtpconfig.get_from_address(), None)
            self.assertEqual(smtpconfig.get_host(), 'localhost')
            self.assertEqual(smtpconfig.get_port(), 25)

        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_config_get_value(self):

            con = configparser.ConfigParser()
            con.add_section(SmtpConfig.DEFAULT)
            smtpconfig = SmtpConfig()

            self.assertEqual(smtpconfig._get_value(con, None, None), None)
            self.assertEqual(smtpconfig._get_value(con, 'blah', None), None)

            self.assertEqual(smtpconfig._get_value(con, SmtpConfig.DEFAULT,
                                                   'yo'), None)

            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_USER,
                    'user')
            self.assertEqual(smtpconfig._get_value(con, SmtpConfig.DEFAULT,
                                                   SmtpConfig.SMTP_USER),
                             'user')

    def test_smtp_config_empty_config_file(self):
        temp_dir = tempfile.mkdtemp()
        try:
            empty_file = os.path.join(temp_dir, 'foo')
            open(empty_file, 'a').close()
            smtpconfig = SmtpConfig(configfile=empty_file)
            self.assertEqual(smtpconfig.get_from_address(), None)
            self.assertEqual(smtpconfig.get_host(), 'localhost')
            self.assertEqual(smtpconfig.get_port(), 25)
        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_config_parse_error(self):
        temp_dir = tempfile.mkdtemp()
        try:
            cfile = os.path.join(temp_dir, 'foo')
            f = open(cfile, 'w')
            f.write('[smtp\nas = xx=12=\n')
            f.flush()
            f.close()
            # this is a coverage test. Basically the
            # _parse_config method called from constructor
            # below should catch the parsing error and
            # not complain
            smtpconfig = SmtpConfig(cfile)
            self.assertEqual(smtpconfig.get_from_address(), None)
            self.assertEqual(smtpconfig.get_host(), 'localhost')
            self.assertEqual(smtpconfig.get_port(), 25)
            self.assertEqual(smtpconfig.get_password(), None)
            self.assertEqual(smtpconfig.get_from_address(), None)
            self.assertEqual(smtpconfig.get_replyto_address(), None)
        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_config_valid_config_file(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(SmtpConfig.DEFAULT)
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_FROM_ADDRESS,
                    'from')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_REPLYTO_ADDRESS,
                    'replyto')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_PASS,
                    'pass')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_PORT,
                    '256')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_HOST,
                    'host')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_USER,
                    'user')
            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()

            smtpconfig = SmtpConfig(configfile=cfile)
            self.assertEqual(smtpconfig.get_user(), 'user')
            self.assertEqual(smtpconfig.get_host(), 'host')
            self.assertEqual(smtpconfig.get_port(), 256)
            self.assertEqual(smtpconfig.get_password(), 'pass')
            self.assertEqual(smtpconfig.get_from_address(), 'from')
            self.assertEqual(smtpconfig.get_replyto_address(), 'replyto')
        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_config_host_and_port_not_set(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(SmtpConfig.DEFAULT)
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_FROM_ADDRESS,
                    'from')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_REPLYTO_ADDRESS,
                    'replyto')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_PASS,
                    'pass')
            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()

            smtpconfig = SmtpConfig(configfile=cfile)
            self.assertEqual(smtpconfig.get_user(), None)
            self.assertEqual(smtpconfig.get_host(), 'localhost')
            self.assertEqual(smtpconfig.get_port(), 25)
            self.assertEqual(smtpconfig.get_password(), 'pass')
            self.assertEqual(smtpconfig.get_from_address(), 'from')
            self.assertEqual(smtpconfig.get_replyto_address(), 'replyto')
        finally:
            shutil.rmtree(temp_dir)

    def test_smtp_config_port_not_a_number(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(SmtpConfig.DEFAULT)
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_PORT,
                    'hello')
            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()

            smtpconfig = SmtpConfig(configfile=cfile)
            self.assertEqual(smtpconfig.get_user(), None)
            self.assertEqual(smtpconfig.get_port(), 25)
        finally:
            shutil.rmtree(temp_dir)

    def test_smtpemailerfactory_nosmtpconfig(self):
        param = D3RParameters()
        fac = SmtpEmailerFactory(param)

        emailer = fac.get_smtp_emailer()
        self.assertEqual(emailer._smtp_host, 'localhost')

    def test_smtpemailerfactory_valid_config(self):
        temp_dir = tempfile.mkdtemp()
        try:

            cfile = os.path.join(temp_dir, 'foo')
            con = configparser.ConfigParser()
            con.add_section(SmtpConfig.DEFAULT)
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_FROM_ADDRESS,
                    'from')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_REPLYTO_ADDRESS,
                    'replyto')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_PASS,
                    'pass')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_PORT,
                    '256')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_HOST,
                    'host')
            con.set(SmtpConfig.DEFAULT, SmtpConfig.SMTP_USER,
                    'user')
            f = open(cfile, 'w')
            con.write(f)
            f.flush()
            f.close()
            param = D3RParameters()
            param.smtpconfig = cfile

            fac = SmtpEmailerFactory(param)

            emailer = fac.get_smtp_emailer()
            self.assertEqual(emailer._user, 'user')
            self.assertEqual(emailer._smtp_host, 'host')
            self.assertEqual(emailer._port, 256)
            self.assertEqual(emailer._password, 'pass')
            self.assertEqual(emailer._fromaddr, 'from')
            self.assertEqual(emailer._replyto, 'replyto')
        finally:
            shutil.rmtree(temp_dir)

    def tearDown(self):
        pass


def try_update_status_from_filesystem(self, task):
    """Runs various tests on update_status_filesystem method
    This is a function cause this needs to be tested on
    all the subclasses as well
    """
    tempDir = tempfile.mkdtemp()
    try:
        task.set_path(tempDir)

        self.assertEqual(task.update_status_from_filesystem(),
                         D3RTask.NOTFOUND_STATUS)

        self.assertEqual(task.get_status(),
                         D3RTask.NOTFOUND_STATUS)

        # Test directory exists but no complete, start, or error files
        task.create_dir()
        self.assertEqual(task.update_status_from_filesystem(),
                         D3RTask.UNKNOWN_STATUS)

        self.assertEqual(task.get_status(),
                         D3RTask.UNKNOWN_STATUS)

        # Test start file exists
        startFile = os.path.join(tempDir, task.get_dir_name(),
                                 'start')
        open(startFile, 'a').close()

        self.assertEqual(task.update_status_from_filesystem(),
                         D3RTask.START_STATUS)

        self.assertEqual(task.get_status(),
                         D3RTask.START_STATUS)

        # Test error file exists
        errorFile = os.path.join(tempDir, task.get_dir_name(),
                                 'error')
        open(errorFile, 'a').close()

        self.assertEqual(task.update_status_from_filesystem(),
                         D3RTask.ERROR_STATUS)

        self.assertEqual(task.get_status(),
                         D3RTask.ERROR_STATUS)

        # Test complete file exists
        completeFile = os.path.join(tempDir, task.get_dir_name(),
                                    'complete')
        open(completeFile, 'a').close()

        self.assertEqual(task.update_status_from_filesystem(),
                         D3RTask.COMPLETE_STATUS)

        self.assertEqual(task.get_status(),
                         D3RTask.COMPLETE_STATUS)

    finally:
        shutil.rmtree(tempDir)


if __name__ == '__main__':
    unittest.main()
