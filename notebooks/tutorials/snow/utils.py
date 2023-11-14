"""
utility function to install API key file
"""
import os


# installation function
def install_api_key(open_topo=False):
    home_dir = os.path.expanduser("~")
    work_dir = os.getcwd()

    # create Topography API key file
    if open_topo:
        topo_key = input("Enter Your OpenTopography API Key: ")
        topo_config_path = os.path.join(work_dir, ".opentopography.txt")

        with open(topo_config_path, "w") as topo_config_file:
            topo_config_file.write(topo_key)
        print(f"OpenTopography API Key file is created at {topo_config_path}.")

    # create CDS API key file
    cds_key = input("Enter Your CDS API Key: ")
    cds_url = "https://cds.climate.copernicus.eu/api/v2"
    cds_config_content = f"url: {cds_url} \nkey: {cds_key}"
    cds_config_path = os.path.join(home_dir, ".cdsapirc")

    with open(cds_config_path, "w") as cds_config_file:
        cds_config_file.write(cds_config_content)

    print(f"CDS API Key file is created at {cds_config_path}")
