from __future__ import print_function

import os
import subprocess
import traceback


print('This is my environment:')
for name, value in os.environ.items():
    print('{name}: {value}'.format(name=name, value=value))

repo_tag = os.environ.get('appveyor_repo_tag', 'false')
tag_name = os.environ.get('appveyor_repo_tag_name', '')
token = os.environ.get('ANACONDA_TOKEN', '')

if repo_tag == 'true' and tag_name.startswith('v'):
    channel = 'main'
else:
    channel = 'dev'
    os.environ['BUILD_STR'] = 'dev'

try:
    resp = subprocess.check_output(['conda', 'build', '--output',
                                    '.conda-recipe'])
except subprocess.CalledProcessError:
    traceback.print_exc()
else:
    file_to_upload = resp.strip().split(os.linesep)[-1]

if not os.path.isfile(file_to_upload):
    print(file_to_upload)
    raise RuntimeError('{name}: not a file'.format(name=file_to_upload))

try:
    subprocess.check_call(['anaconda', '-t', token, 'upload', '--force',
                           '--user', 'landlab', '--channel', channel,
                           file_to_upload])
except subprocess.CalledProcessError:
    traceback.print_exc()
