import numpy as np
import shutil
from astroquery.alma.core import Alma as ALMA
import sys 
version = sys.version_info[0]
if  version == 3:
    import pathlib
elif version == 2:
    import pathlib2 as pathlib
    
class AlmaQuery:
    def __init__(self, dict):
        self.query = Alma.query(dict)
        
        self.PI = self.query['PI name'][0].split(',')[0]
        self.project_filelinks = Alma.stage_data(self.query[0]['Member ous id'])
        self.download_link = self.project_filelinks[self.project_filelinks['size'].argmin()]
    
    def filter_links(self, string):
        for link in self.project_filelinks:
            if string in link['URL']:
                self.download_link = link
        print('\nSize = {} GB'.format(self.download_link['size']))
            
    def url(self):
        try:
            url = self.download_link['URL'].replace('dataPortal', 'rh')
        except AttributeError:
            url = self.project_filelinks[0]['URL'].replace('dataPortal', 'rh')
        url = str(url[:url.find('ALMA')])
        print(url)
    
    def download(self):
        # download, extract README and fits files, and then delete everything else
        filepaths = Alma.download_and_extract_files(self.download_link['URL'], regex=r'.*\.fits$|.*README$')            
        
        # move downloaded files to the ALMA_data directory
        self.name = self.query['Source name'][0]
        self.directory = pathlib.Path('ALMA_data/' + self.name + '_' + self.PI)
        self.directory.mkdir(exist_ok=True)
        for path in filepaths:
            current_path = pathlib.Path(path)
            new_path = self.directory.joinpath(current_path.name)
            shutil.move(str(current_path), str(new_path))
                
Alma = ALMA()
Alma.login('cdaley')
Alma.cache_location = 'ALMA_cache'


# q = AlmaQuery({'source_name_resolver' : 'TW Hya', 'project_code' :'2016.1.00629.S'})
# q.url()
# if version == 2:
#     q.download()

q = AlmaQuery({'source_name_alma' : 'Elia_2-27', 'project_code' :'2013.1.00498.S'})
if version == 2:
    q.download()
    
# for link in q.download_links:
#     if not all([string in link['URL'] for string in ['2016.1.01042.S_uid___A001_X892_X49_001_of_001.tar']]):
#         print(link)
# q = AlmaQuery({'project_code' : '2016.1.00195.S'})
# if version == 2:
#     q.download()
