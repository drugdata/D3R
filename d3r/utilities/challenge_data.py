__author__ = 'j-wags'

import logging
import glob
import os

class ChallengeData(object):
    """
    A class for interacting with challenge data packages
    """
    
    def __init__(self, get_challenge_data_dir):
        self.get_challenge_data_dir = get_challenge_data_dir
        self.abs_get_challenge_data_dir = os.path.abspath(get_challenge_data_dir)
        
        self.weekly_challenge_packages = glob.glob('%s/celpp_week*/' %(self.abs_get_challenge_data_dir))
        self.weekly_challenge_folder_names = []
        print 'self.weekly_challenge_packages',self.weekly_challenge_packages

        self.week_challenge_dict = {}
        for wcp in self.weekly_challenge_packages:
            print 'wcp', wcp
            wcp_folder_name = wcp.strip('/').split('/')[-1]
            self.weekly_challenge_folder_names.append(wcp_folder_name)
            print 'wcp_folder_name', wcp_folder_name
            targets = glob.glob('%s/????/' %(wcp))
            print 'targets', targets
            self.week_challenge_dict[wcp_folder_name] = targets
            
        

    def get_week_names(self):
        """
        Gets the bottom-level folder names in this challenge pacakge ("celpp_weekXX_XXXX"), with slashes removed
        :return: list of strings
        """
        return self.weekly_challenge_folder_names
    
    def get_weekly_challenge_packages(self):
        """
        A function to return the all challenge data packages inside of this directory. Expects subdirectory contining a folder named, for example, celpp_week25_2016. Returns a list of absolute folder paths
        :return: A list of strings, each an absolute path to a challenge data package
        """
        
        return self.weekly_challenge_packages
                         


    def is_valid_for_celpp(self):
        """
        Determines whether the challenge data package contains exactly one week with at least one target
        :return: boolean
        """
        
        # Ensure that exactly one week is in here
        if len(self.week_challenge_dict.keys()) != 1:
            logging.info("Challenge package contains more than 1 week and is therefore not suitable for CELPP (contains %r)" %(self.week_challenge_dict.keys()))
            return False

        # Ensure that the week has at least one target
        if len(self.week_challenge_dict[self.week_challenge_dict.keys()[0]]) == 0:
            logging.info("Week in challenge package does not contain any targets (week folder %s)" %(self.week_challenge_dict[self.week_challenge_dict.keys()[0]]))
            return False

        return True

    def get_targets(self):
        """
        Returns the {week:[targets]} dictionary for this challenge package. If is_valid_for_celpp is True, this dictionary will have one key, with the corresponding value being a list of strings of absolute directories of targets.
        :return: dict
        """
        
        return self.week_challenge_dict
        
