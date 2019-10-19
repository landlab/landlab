#! /usr/bin/env python


import os
import re
import subprocess
import sys
from collections import OrderedDict, defaultdict

import jinja2

CHANGELOG = """
# Change Log
All notable changes to landlab will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](https://semver.org/).

This file was auto-generated using `scripts/make_changelog.py`.

{% for tag, sections in releases.iteritems() %}
## [{{ tag }}] {{ release_date[tag] }}
{% for section, changes in sections.iteritems() %}
### {{section}}
{% for change in changes -%}
* {{ change }}
{% endfor -%}
{% endfor -%}
{% endfor -%}
""".strip()

SECTIONS = ['Added', 'Changed', 'Deprecated', 'Removed', 'Fixed', 'Security']


def git_log(start=None, stop='HEAD'):
    cmd = ['git', 'log', '--first-parent', 'master', '--merges',
           '--topo-order',
           # '--pretty=message: %s+author:%an+body: %b'],
           '--pretty=%s [%an]',
           # '--oneline',
          ]
    if start:
        cmd.append('{start}..{stop}'.format(start=start, stop=stop))
    return subprocess.check_output(cmd).strip()


def git_tag():
    return subprocess.check_output(['git', 'tag']).strip()


def git_tag_date(tag):
    return subprocess.check_output(['git', 'show', tag,
                                    '--pretty=%ci']).strip().split()[0]


def releases(ascending=True):
    tags = git_tag().split(os.linesep) + ['HEAD']
    if ascending:
        return tags
    else:
        return tags[::-1]


def format_pr_message(message):
    m = re.match(
        'Merge pull request (?P<pr>#[0-9]+) '
        'from (?P<branch>[\S]*)'
        '(?P<postscript>[\s\S]*$)',
        message
    )
    if m:
        return '{branch} [{pr}]{postscript}'.format(**m.groupdict())
    else:
        raise ValueError('not a pull request')


def format_changelog_message(message):
    m = re.match(
        '(?P<first>\w+)(?P<rest>[\s\S]*)$',
        message
    )
    word = m.groupdict()['first']
    if word in ('Add', 'Fix', 'Deprecate', ):
        return word + 'ed' + m.groupdict()['rest']
    elif word in ('Change', 'Remove', ):
        return word + 'd' + m.groupdict()['rest']
    else:
        return message


def prettify_message(message):
    if message.startswith('Merge branch'):
        return None

    try:
        message = format_pr_message(message)
    except ValueError:
        message = format_changelog_message(message)
    return message


def brief(start=None, stop='HEAD'):
    changes = []
    for change in git_log(start=start, stop=stop).split(os.linesep):
        if change:
            message = prettify_message(change)
            if message:
                changes.append(message)

    return changes


def group_changes(changes):
    groups = defaultdict(list)
    for change in changes:
        if change.startswith('Add'):
            group = 'Added'
        elif change.startswith('Deprecate'):
            group = 'Deprecated'
        elif change.startswith('Remove'):
            group = 'Removed'
        elif change.startswith('Fix'):
            group = 'Fixed'
        elif change.startswith('Security'):
            group = 'Security'
        else:
            group = 'Changed'
        groups[group].append(change)
    return groups


def main():
    tags = releases(ascending=False)
    changelog = OrderedDict()
    release_date = dict()
    for start, stop in zip(tags[1:], tags[:-1]):
        changes = brief(start=start, stop=stop)
        if changes:
            changelog[stop] = group_changes(changes)
            release_date[stop] = git_tag_date(stop)

    env = jinja2.Environment(loader=jinja2.DictLoader({'changelog': CHANGELOG}))
    print(env.get_template('changelog').render(releases=changelog,
                                               release_date=release_date))


if __name__ == '__main__':
    sys.exit(main())
