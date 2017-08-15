"""
--------------------------------------------------------------------------------
FILE: Application.py
AUTH: thorsilver
UPDT: digimus
VERS: 1.1 July 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""

from __future__ import print_function # ensure Python 3.x print form known
from math import tanh

# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def _calculate_grant_quality(author, params, rng):
    """
    Calculate quality of grant.

    Three components:
    - time invested in grant writing
    - research output to date
    - (optionally) applicants research quality
    - random noise
    """
    if author.time_grant <= 0.0:
        print("Postdoc?")
    quality = params['weight_grant'] * \
              tanh(params['grant_slope'] * author.time_grant) + \
              (1 - params['weight_grant']) * \
              tanh(params['research_slope'] * author.research_sum)
    if params['rq_counts']:
        quality += author.research_quality
    noise = rng.normalvariate(0.0, params['grant_noise'])
    quality *= (1.0 + noise)
    return quality

# --------------------------------------------------
# Main class
# --------------------------------------------------

class Application(object):

    """
    A grant application.

    Attributes:
    - author_id : the id of the grant author
    - author_time : amount of time invested by the author
    - author_quality : research_quality of authoring academic
    - grant_quality : the quality of the grant, used for ranking
    """

    def __init__(self, author, params, rng):
        self.author_id = author.ident
        self.author_time = author.time_grant
        self.author_quality = author.research_quality
        self.grant_quality = _calculate_grant_quality(author, params, rng)

    def __cmp__(self, other):
        """
        Simulate the deprecated cmp function to return -1 if a<b, 0 if a==b,
        or 1 otherwise. This replaces the original call of
        cmp(other.grant_quality, self.grant_quality)
        """

        if other.grant_quality < self.grant_quality:
            return -1
        elif other.grant_quality == self.grant_quality:
            return 0
        return 1

    def dummy1(self):

        """
        *** Dummy method to ensure there is a method for this class ***
        *** Remove this later ***
        """
        return self.grant_quality

    def dummy2(self):

        """
        *** Dummy method to ensure there is a method for this class ***
        *** Remove this later ***
        """
        return self.grant_quality


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
