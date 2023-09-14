"""
utility function to install API key file
"""
import os


# installation function
def install_api_key():
    home_dir = os.path.expanduser('~')

    # create CDS API key file
    cds_key = input('Enter Your CDS API Key: ')
    cds_url = 'https://cds.climate.copernicus.eu/api/v2'
    cds_config_content = 'url: {} \nkey: {}'.format(cds_url, cds_key)
    cds_config_path = os.path.join(home_dir, '.cdsapirc')

    with open(cds_config_path, 'w') as cds_config_file:
        cds_config_file.write(cds_config_content)

    print('CDS API Key file is created at {}'.format(cds_config_path))
