"""
Docstring
"""

import json
import os


def _init_parms():
    """
    Initialise parameters at start of run. Values are taken from
    ./_default.json file
    """

    with open('json/_default.json') as json_data:
        parms = json.load(json_data)

    return parms

def _init_prefs():
    """
    Initialise preferences from the .majs file in users' home directory
    """

    pref_file = os.path.expanduser('~/.majs')
    with open(pref_file) as json_data:
        prefs = json.load(json_data)

    return prefs

class Dummy(object):
    """
    Docstring
    """

    def __init__(self):
        """
        Initialise parameters and preferences
        """

        self.parms = _init_parms()
        self.prefs = _init_prefs()


    def testing(self):
        """
        Try out
        """

        print(parms)
        print('\n\n')
        print(prefs)
