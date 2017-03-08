# -*- coding: utf-8 -*-

__author__ = 'churas'

import logging
import re

logger = logging.getLogger(__name__)


class Participant(object):
    """Represents a an external participant in the challenge
    """

    def __init__(self, name, d3rusername, guid, email):
        """Constructor
        """
        self._name = name
        self._d3rusername = d3rusername
        self._guid = guid
        self._email = email

    def get_name(self):
        """Gets name
        """
        return self._name

    def get_d3rusername(self):
        """Gets d3rusername
        """
        return self._d3rusername

    def get_guid(self):
        """Gets guid
        """
        return self._guid

    def get_email(self):
        """Gets email
        """
        return self._email


class ParticipantDatabase(object):
    """Database object that contains a set of `Participants`
    """

    def __init__(self, participants):
        """Constructor
        :param participants: List of participants
        """
        self._participants = participants

    def get_participant_by_guid(self, guid,
                                exact_match=False):
        """Searches database for `Participant` with matching `guid`
        :param guid: guid or unique id for `Participant`
        :param exact_match: if set to True only guid matching
                            exactly with a `Participant` will be returned.
                            If False then the method will try for an exact
                            match, but if that fails the first `Participant`
                            that matches guid if _# is tripped of is used.
        :returns `Participant` with matching `guid` or None if not found
        """
        if guid is None:
            logger.warning('guid passed in was None')
            return None

        if self._participants is None:
            logger.warning('participants is None')
            return None

        for p in self._participants:
            if p.get_guid() == guid:
                return p

        if exact_match is True:
            logger.debug('No match found and exact_match is set '
                         'to True. returning None for guid: ' + guid)
            return None

        # stripped_guid = re.sub('_[0-9A-Za-z-]+$', '', guid)
        possible_guids = re.findall('([0-9]{5})_[0-9A-Za-z-]+$', guid)
        if len(possible_guids) > 1:
            logger.debug('Stripped guid is not unique. Stripping ' + guid +
                         ' yielded ' + possible_guids)
            return None
        if len(possible_guids) == 0:
            logger.debug('guid stripping did not yield a match. Original'
                         ' guid was ' + guid)
            return None
        stripped_guid = possible_guids[0]
        if stripped_guid != guid:
            logger.debug('Stripped guid is different: ' + stripped_guid +
                         'performing search again')
            for p in self._participants:
                if p.get_guid() == stripped_guid:
                    return p

        return None

    def get_participants(self):
        """Gets participants as list
        :returns: List of `Participants`
        """
        return self._participants


class ParticipantDatabaseFromCSVFactory(object):
    """Factory class that generates `ParticipantDatabase` object from CSV file
    """

    def __init__(self, csvfile):
        """Constructor that takes CSV file `csvfile` to generate
        `ParticipantDatabase` object.  The format of the CSV file
        is as follows:

        name,d3rusername,guid,email
        bob smith,bsmith,8675309,bob@bob.com

        :param csvfile: Path to `Participant` CSV file
        """
        self._csvfile = csvfile

    def get_participant_database(self):
        """Parses CSV file set in constructor to create `ParticipantDatabase`
        :returns `ParticipantDatabase` from CSV file or None if there was an
        error
        """
        if self._csvfile is None:
            logger.warning('No csv file set')
            return None
        try:
            f = open(self._csvfile, 'rU')
            counter = 0
            plist = []
            for line in f:
                if counter is 0:
                    if line.startswith('name'):
                        logger.debug('Skipping header line')
                        counter = + 1
                        continue

                splitline = line.rstrip().split(',')
                if len(splitline) is not 4:
                    logger.warning('Problems splitting line ' + line +
                                   ' got ' + str(len(splitline)) +
                                   ' elements expecting 4')
                    counter = + 1
                    continue
                plist.append(Participant(splitline[0].strip(),
                                         splitline[1].strip(),
                                         splitline[2].strip(),
                                         splitline[3].strip()))
                counter = + 1
            f.close()
            return ParticipantDatabase(plist)
        except Exception:
            logger.exception('Caught exception')
        return None
