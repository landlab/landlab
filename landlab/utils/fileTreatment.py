import os


def check_out_dirs(out_file_dir):

    if os.path.isdir(out_file_dir) is False:
        os.makedirs(out_file_dir)
        print("New Directory created for output files:")
        print(out_file_dir)