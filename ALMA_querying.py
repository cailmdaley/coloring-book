import numpy as np
import pathlib2
import shutil
from astroquery.alma.core import Alma
Alma.cache_location = 'ALMA_cache'

class AlmaQuery:
    def __init__(self, name, code=None, download=False):
        self.query = Alma.query_object(name)
        if code is not None:
            self.project = self.query[self.query['Project code'] == '2015.1.00686.S']
            self.name = self.project['Source name'][0]
        if download is True:
            # choose smallest of the tar files to download (should be data product)
            filelink_list = Alma.stage_data(self.project['Member ous id'])
            filelink = filelink_list[filelink_list['size'].argmin()]
            print('\nSize = {} GB'.format(filelink['size']))
            
            # download, extract README and fits files, and then delete everything else
            filepaths = Alma.download_and_extract_files(filelink['URL'], regex=r'.*\.fits$|.*README$')            
            
            # move downloaded files to the ALMA_data directory
            self.directory = pathlib2.Path('ALMA_data/' + self.name)
            self.directory.mkdir(exist_ok=True)
            for path in filepaths:
                current_path = pathlib2.Path(path)
                new_path = self.directory.joinpath(current_path.name)
                shutil.move(str(current_path), str(new_path))
                
            
q = AlmaQuery('TW Hya', '2015.1.00686.S', download=True)
