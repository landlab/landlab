from __future__ import print_function

import os
import sys
import subprocess
import traceback


# print('This is my environment:')
# for name, value in os.environ.items():
#     print('{name}: {value}'.format(name=name, value=value))

print('Using python: {prefix}'.format(prefix=sys.prefix))

repo_tag = os.environ.get('APPVEYOR_REPO_TAG', 'false')
tag_name = os.environ.get('APPVEYOR_REPO_TAG_NAME', '')
token = os.environ.get('ANACONDA_TOKEN', 'NOT_A_TOKEN')

if repo_tag == 'true' and tag_name.startswith('v'):
    channel = 'main'
else:
    channel = 'dev'
    os.environ['BUILD_STR'] = 'dev'

print('Uploading to {channel} channel'.format(channel=channel))

try:
    cmd = ' '.join(['conda', 'build', '--output', '.conda-recipe'])
    resp = subprocess.check_output(cmd, shell=True)
except subprocess.CalledProcessError:
    traceback.print_exc()
else:
    file_to_upload = resp.strip().split(os.linesep)[-1]

print(file_to_upload)

if not os.path.isfile(file_to_upload):
    raise RuntimeError('{name}: not a file'.format(name=file_to_upload))

try:
    cmd = ' '.join(['anaconda', '-t', token, 'upload', '--force', '--user',
                    'landlab', '--channel', channel, file_to_upload])
    subprocess.check_call(cmd, shell=True)
except subprocess.CalledProcessError:
    traceback.print_exc()
