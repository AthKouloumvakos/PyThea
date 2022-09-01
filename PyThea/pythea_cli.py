'''
When PyThea is installed in an enviroment,
with this script you can directry call the package as an excecutable.
$ python3 -m pip install ./ --use-feature=in-tree-build
$ pythea streamlit
'''
import os
import sys
from typing import Optional

import click
import streamlit
import streamlit.web.bootstrap as bootstrap
from streamlit.runtime.credentials import check_credentials

from .version import version as _vesrsion


@click.group()
def main():
    pass


@main.command('streamlit')
def main_streamlit():
    """Run the main program in browser."""

    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'PyThea_app.py')

    _main_run(filename, flag_options=None)


def _main_run(file, args=None, flag_options=None):
    if args is None:
        args = []

    if flag_options is None:
        flag_options = {}

    command_line = _get_command_line_as_string()

    # Set a global flag indicating that we're "within" streamlit.
    streamlit._is_running_with_streamlit = True

    check_credentials()

    bootstrap.run(file, command_line, args, flag_options)


def _get_command_line_as_string() -> Optional[str]:
    import subprocess

    parent = click.get_current_context().parent
    if parent is None:
        return None

    if 'streamlit.cli' in parent.command_path:
        raise RuntimeError(
            'Running streamlit via `python -m streamlit.cli <command>` is'
            ' unsupported. Please use `python -m streamlit <command>` instead.'
        )

    cmd_line_as_list = [parent.command_path]
    cmd_line_as_list.extend(sys.argv[1:])
    return subprocess.list2cmdline(cmd_line_as_list)


@main.command('help')
@click.pass_context
def help(ctx):
    """Print this help message."""
    # For 'pythea --help' instead of 'pythea help'.
    import sys
    sys.argv[1] = '--help'
    main()


@main.command('version')
def version():
    """Print PyThea's version number."""
    print(f'PyThea installed version is: {_vesrsion}')


@main.command('docs')
def docs():
    """Show doc page in browser."""
    print('Showing help page in browser...')
    from streamlit import util

    util.open_browser('https://www.pythea.org/en/docs/')


@main.command('update')
def update():
    """Update PyThea."""
    print('Not implemented yet.')

# TODO Include the test here.


if __name__ == '__main__':
    main()
