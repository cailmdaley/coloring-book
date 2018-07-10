from astroquery.alma.core import Alma
import numpy as np
import logging 
logging.basicConfig(format='%(message)s', level=logging.WARNING)

class AlmaQuery:
    def __init__(self, name, code=None, download=False):
        self.query = Alma.query_object(name)
        if code is not None:
            self.project = self.query[self.query['Project code'] == '2015.1.00686.S']
        if download is True:
            filelinke_list = Alma.stage_data(self.project['Member ous id'])
            filelink = filelinke_list[filelinke_list['size'].argmin()]
            print('\nSize = {} GB'.format(filelink['size']))

q = AlmaQuery('TW Hya', '2015.1.00686.S', download=True)
