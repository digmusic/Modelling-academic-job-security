"""
--------------------------------------------------------------------------------
FILE: FundingAgency.py
AUTH: thorsilver
VERS: 1.0 March 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""

class FundingAgency(object):

    """
    The funding body.

    Collects applications and allocates grants.  Typically there will only
    be one instance of this per simulation.

    Agents are sorted into a given number of pools, and the top agents from
    each pool are awarded grants.  If the number of pools is one, then agents
    are competing against the entire population.
    """

    def __init__(self, params):
        self.pools = []
	# allocate individuals to pools
        for _ in range(params['grant_pools']): # *** was xrange ***
            self.pools.append([])
	# store data on successful applications
        self.successful_app_stats = []
        self.num_grants = params['starting_grant_fund']


    def _get_grant_recipients(self, params, size, scale, limited):

        """
        Helper to return a list of authors whose applications were successful.
        """

        if limited:
            num_grants = self.num_grants
        else:
            num_grants = params['grant_proportion'] * scale * size
            print('Number of grants available: {}').format(num_grants)

        per_pool = [int(num_grants / params['grant_pools'])] * \
                params['grant_pools']

        leftover = int(num_grants % params['grant_pools'])

        for item in range(leftover):
            per_pool[item] += 1

        success = []
        for item in range(params['grant_pools']):
            success.extend([app for app in self.pools[item][:per_pool[item]]])
            self.successful_app_stats.extend([(
                app.author_quality,
                app.author_time) for app in self.pools[item][:per_pool[item]]])

        return [app.author_id for app in success]


    def add_application(self, new_app, rng):

        """
        Add new_app to a randomly chosen pool.
        """

        self.pools[rng.randint(0, len(self.pools)-1)].append(new_app)


    def clear_applications(self):

        """
        Clear all current applications.
        """

        for item in range(len(self.pools)):
            self.pools[item] = []


    def rank_applications(self):

        """
        *** Description ***
        """

        # *** Need to specify a key to sort on? Currently complains of
        # TypeError: '<' not supported between instances of
        # 'Application' and 'Application' ***

        for pool in self.pools:
            pool.sort(key=lambda Application: Application.grant_quality)

    def get_grant_recipients(self, params, size):

        """
        Return a list of authors whose applications were successful.
        """
        return self._get_grant_recipients(params, size, 0.75, False)


    def get_grant_recipients_pdr(self, params, size):

        """
        Return a list of authors whose applications were successful.
        """
        return self._get_grant_recipients(params, size, 1.0, False)


    def init_grants(self, params):

        """
        *** Description ***
        """
        self.num_grants = params['starting_grant_fund']


    def update_grants(self, params):

        """
        *** Description ***
        """
        self.num_grants += self.num_grants * params['yearly_increase']
        self.num_grants = int(self.num_grants)


    def get_recipients_limited(self, params, size):

        """
        Return a list of successful applications in limited funding case.
        """
        return self._get_grant_recipients(params, size, 1.0, True)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
