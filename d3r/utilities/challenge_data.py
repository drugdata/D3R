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
        
        self.weekly_challenge_packages = glob.glob('%s/celpp_week*' %(self.abs_get_challenge_data_dir))
        
        self.week_challenge_dict = {}
        for wcp in self.weekly_challenge_packages:
            wcp_folder_name = wcp.strip('/').split('/')[-1]
            targets = glob.glob('%s/????/')
            self.week_challenge_data[wcp_folder_name] = targets
            
        

    
        
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
        if len(self.week_challenge_data.keys()) != 1:
            return False

        # Ensure that the week has at least one target
        if len(self.week_challenge_data[self.week_challenge_data.keys()[0]]) == 0:
            return False

        return True

    def get_targets(self):
        """
        Returns the {week:[targets]} dictionary for this challenge package. If is_valid_for_celpp is True, this dictionary will have one key, with the corresponding value being a list of strings of absolute directories of targets.
        :return: dict
        """
        
        return self.week_challenge_data
        
