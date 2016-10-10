__author__ = 'j5wagner'

import re

class ReadText(object):
    """
    CELPP <pdbid>.txt file reader
    """
    ## A dictionary of regular expressions for parsing target txt file lines
    PARSING_DICT = {'query':         '(\S+)',
                    'ph':            '(\S+)',
                    'ligand':        '([A-Za-z0-9]+)',
                    'inchi':         '(InChI=.*)',
                    'size':          '([0-9]+)',
                    'rotatable_bond':'([0-9]+)',
                    'LMCSS':         (['cand_id','lig_name','chain','size','mcss_size','resolution'],
                                      '([a-zA-Z0-9]{4}), ([a-zA-Z0-9]+), chain: (.), [(]size:[ ]([0-9]+), mcss_size:[ ]([0-9]+), resolution:[ ]*([0-9.]+)[)]\w*'),

                    'SMCSS':         (['cand_id','lig_name','chain','size','mcss_size','resolution'],
                                      '([a-zA-Z0-9]{4}), ([a-zA-Z0-9]+), chain: (.), [(]size:[ ]*([0-9]+), mcss_size:[ ]*([0-9]+), resolution:[ ]*([0-9.]+)[)]\w*'),

                    'hiResHolo':     (['cand_id','lig_name','chain','resolution'],
                                      '([a-zA-Z0-9]{4}), ([a-zA-Z0-9]+), chain: (.), [(]resolution:[ ]*([0-9.]+)[)]\w*'),

                    'hiTanimoto':    (['cand_id','lig_name','chain','tanimoto_similarity', 'resolution'],
                                      '([a-zA-Z0-9]{4}), ([a-zA-Z0-9]+), chain: (.), [(]tanimoto_similarity:[ ]*([0-9.]+), resolution:[ ]*([0-9.]+)[)]\w*'),

                    'hiResApo':      (['cand_id'],
                                      '([a-zA-Z0-9]{4})\w*'),
                    }
    
    ## A list of line types that return a dict instead of a single value
    DICT_RETURN_KEYS = ['LMCSS', 
                        'SMCSS',
                        'hiResHolo',
                        'hiTanimoto',
                        'hiResApo']

    def parse_line(self, line):
        """ 
        Parses a single line of a CELPP target .txt file
        :param: the line for processing
        :return: string - the key for the given line
        :return: string or dict - the value for the given line
        """
        line_split = line.split(',')
        txt_key = line_split[0]
        if txt_key[0] == '#':
            #print 'Skipping commmented line %s' %(line)
            return None, None

        if not(txt_key in self.PARSING_DICT.keys()):
            #raise Exception('Unable to parse line %r: key %s is invalid' %(line, key))
            #print 'Unable to parse line %r: key %s is invalid' %(line, txt_key)
            raise Exception('Unable to parse line %r: key %s is invalid' %(line, txt_key))
            return None, None

        value = ','.join(line_split[1:])
        #print
        #print 
        #print 'line', line
        if txt_key in self.DICT_RETURN_KEYS:
            data_keys = self.PARSING_DICT[txt_key][0]
            parsing_re = self.PARSING_DICT[txt_key][1]

            #print 'data_keys', data_keys
            #print 'parsing_re', parsing_re
            parsed_value = re.findall(parsing_re, value)
            #print 'parsed_value', parsed_value

            # Fix layering problem for candidates with only one field
            if not(type(parsed_value[0]) is tuple):
                parsed_value = tuple([parsed_value])
    

            if len(parsed_value[0]) != len(data_keys):
                #print "Didn't parse as many data values as data keys"
                raise Exception("Didn't parse as many data values as data keys")
                return None, None
            return_dict = {}
            for data_key, value in zip(data_keys,parsed_value[0]):
                #if data_key == 'other_info':
                    # Parse the value. 
                return_dict[data_key] = return_dict.get(data_key,[]) + [value]
            return txt_key, return_dict

        else:
            parsing_re = self.PARSING_DICT[txt_key]
            #print 'parsing_re', parsing_re
            parsed_value = re.findall(parsing_re, value)
            if len(parsed_value) != 1:
                #print "Parsing %s using RE %s returned %r. This parsing must return one unique match. Skipping line." %(value, parsing_re, parsed_value)
                raise Exception("Parsing %s using RE %s returned %r. This parsing must return one unique match. Skipping line." %(value, parsing_re, parsed_value))
                return None, None
            return txt_key, parsed_value[0]

    def parse_txt(self, txt_name):
        """Parses a CELPP target .txt file, returning a
        dictionary. Normal lines (eg "ph") are parsed into a
        key:[value] pair (eg. 'ph':['7.4']). As candidate information
        lines have a variable amount of data, they are parsed into a
        dictionary (eg. 'hiTanimoto':[{'cand_id': ['4di5'],
        'lig_name': ['DPO'], 'resolution': ['2.3'],
        'tanimoto_similarity': ['1.0'], 'chain': ['A']}] 
        :param: txt_name - the name of the txt file to be parsed 
        :return: a dictionary of all the information from the txt file.
        """

        data_dict = {}
        with open(txt_name) as fo:
            for line in fo:
                key, parsed_value = self.parse_line(line)
                if key is None:
                    continue
                    
                data_dict[key] = data_dict.get(key,[]) + [parsed_value]

        #print 
        #print
        #print
        #for key in data_dict.keys():
            #print key, ':', data_dict[key]
        return data_dict
                    

