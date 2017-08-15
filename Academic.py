"""
--------------------------------------------------------------------------------
FILE: Academic.py
AUTH: thorsilver
UPDT: digimus
VERS: 1.1 July 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""
from random import randint

# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def _get_update_step(params, rng):
    """
    Get update step size.
    """

    if params['self_update_width_fixed']:
        return params['self_update_width']
    return abs(rng.gauss(0.0, params['self_update_width']))

def _res1(parm, quality):
    """
    Calc_research helper function (type 1)
    """
    return (1.0 - parm) * quality

def _res2(parm, time, held, bonus, quality):
    """
    Calc_research helper function (type 2)
    """
    return ((1.0 + parm - time) + (held * bonus)) * quality

def _res3(time, held, bonus, quality):
    """
    Calc_research helper function (type 3)
    """
    return ((1.0 - time) + (held * bonus)) * quality

# --------------------------------------------------
# Main class
# --------------------------------------------------

class Academic(object):

    """
    An academic agent.

    Attributes:
    - ident : unique identifier
    - research_quality : quality of research (random float in range 0 - 1)
    - applying: True if agent is currently applying for grants
    - grant_held : grant held for current year (boolean 1/0)
    - grant_count : number of grants received over total life
    - time_grant : time allocated to grant writing (initially random)
    - research : research produced this year
    - research_sum : cumulative research output (initially 0)
    - memory : memory, used for strategy update rules
    - postdoc_status : is the agent a postdoc (boolean 1/0)
    - former_postdoc : set to 1 if agent became permanent (boolean 1/0)
    - exited_postdoc : set to 1 if agent can't find job and exits (boolean 1/0)
    - num_postdocs : number of contracts over lifetime
    - contract_length : contract length remaining in semesters (max 10)
    - newb : are they a postdoc newb? (0.4 reduction in research time)

    Time allocated to research is currently assumed to be 1.0 - time_grant
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self, ident, params, rng):

        self.applying = True
        self.career_length = params['career_end']
        self.contract_length = 100
        self.exited_postdoc = 0
        self.former_postdoc = 0
        self.grant_count = 0    # total number of grants held
        self.grant_held = False
        self.ident = ident      # renamed by IW
        self.made_redundant = False
        self.memory = []        # grows at tail
        self.newb = 0
        self.num_postdocs = 0
        self.params = params
        self.postdoc_status = 0
        self.research = 0.0
        self.research_quality = rng.random()
        self.research_sum = 0.0
        self.retired = False
        self.set_random_time_grant(params['init_time'], rng)
        self.tenured = True
        self.time_grant = 0     # added by IW



    # --------------------------------------------------
    # Independent methods
    # --------------------------------------------------

    def check_time_bounds(self):
        """
        Ensure that an agents time_grant allocation is between 0 and 1.
        """

        if self.time_grant > 1.0:
            self.time_grant = 1.0
        elif self.time_grant < 0.0:
            self.time_grant = 0.0

        # must make a minimum effort to apply; time_grant can only be 0 if
	# applying is False.
        if self.applying and self.time_grant < 0.1:
            self.time_grant = 0.1


    def get_mean_research(self, iters):

        """
	Get agent's mean research output over the preceding t iterations.

	(or over less than t, if fewer iterations have occurred)
	"""

        output_values = [m[2] for m in self.memory[-iters:]]
        if output_values == []:
            return 0.0
        return sum(output_values) / len(output_values)


    def retire(self):
        """
        Check the career length of the agent. If they have been in post
        for 10 years or more there's a 20% chance in a given year they
        may choose to retire or leave the profession.

        """
        if self.career_length >= 10 and randint(1, 100) <= 20:
            self.made_redundant = True
            self.retired = True
            self.applying = False

    def set_random_time_grant(self, time_range, rng):
        """
        Randomise the time allocated to grant writing.
        Current modification limits to 0.1 increments.
        """
        # if self.postdoc_status == 1:
        #    self.time_grant = 0.0
        # else:
        self.time_grant = rng.random() / (1.0 / time_range)
        self.time_grant = int(self.time_grant * 10) / 10.0
        # agent must invest a minimum of 10% of their time
        # (unless they drop out)
        if self.time_grant < 0.1:
            self.time_grant = 0.1


    def update_memory(self, params):
        """
        Update an agent's memory with current strategy, success and output.
        """

        self.memory.append((self.time_grant, self.grant_held, self.research))
        if len(self.memory) > params['memory_size']:
            del self.memory[0]


    # --------------------------------------------------
    # Level 0 methods: dependent upon helpers
    # --------------------------------------------------

    def check_reentry(self, params, rng):
        """
        Possibly reintroduce a previously dropped-out agent.
        """

        if rng.random() < params['prob_reentry']:
            self.set_random_time_grant(params['reentry_range'], rng)
            self.applying = True
            self.memory = []

    def calc_research(self, time_grant, grant_held, grant_bonus,
                      research_quality):
        """
        Calculate agent research quality.
        """

        if self.params['growing_pop'] == 1:
            if self.newb >= 1:
                return _res1(self.params['newb_time'], research_quality)
            elif self.grant_held:
                return _res2(self.params['manager_penalty'], time_grant,
                             grant_held, grant_bonus, research_quality)
        elif self.params['use_postdocs'] == 1:
            if self.newb >= 1:
                return _res1(self.params['newb_time'], research_quality)
            elif self.contract_length <= 2 and self.postdoc_status == 1:
                return _res1(self.params['jobhunt_time'], research_quality)
            elif self.grant_held:
                return _res2(self.params['manager_penalty'], time_grant,
                             grant_held, grant_bonus, research_quality)
            elif self.former_postdoc == 1 and self.params['mentored_pdrs'] == 1:
                return _res2(self.params['mentoring_bonus'], time_grant,
                             grant_held, grant_bonus, research_quality)
        # for anything else use _res3
        return _res3(time_grant, grant_held, grant_bonus, research_quality)

    # --------------------------------------------------
    # Level 1 methods
    # --------------------------------------------------

    def produce_research(self, params):

        """
	Calculate an agent's research output for a single year.
	Also updates research_sum and memory.
	"""

        self.research = self.calc_research(
            self.time_grant, self.grant_held, params['grant_bonus'],
            self.research_quality)
        self.research_sum += self.research
        self.update_memory(params)
        return self.research

    ### SELF LEARNING RULES ###

    def update_strategy_self_memory(self, params, rng):

        """
        Update agent strategy based on self learning.

        1. agents who do not hold a grant, but who have held a grant recently
            increase their time allocation
        2. agents who do not hold a grant, and who have not held a grant
            recently stop applying (invest no time)
        3. agents who hold a grant and have consistently held a grant in
            recent iterations decrease their time allocation
        4. agents who hold a grant, but who have not held a grant recently
            hold their time constant

        Agents who have previously stopped applying have a probability of
        re-entering competition at each iteration.

        WCSS paper - MEMORY model (with various memory / re-entry parameters)
        """

        step = _get_update_step(params, rng)
        recent = [x[1] for x in self.memory]
        runl = params['run_length']

        if not self.grant_held:
            if len(recent) > runl and recent[-runl:] == [False]*runl:
                self.time_grant = 0.0
                self.applying = False
            else:
                #self.time_grant *= (1.0 + step)
                self.time_grant += step
        else:
            if len(recent) > runl and recent[-runl:] == [True]*runl:
                #self.time_grant *= (1.0 - step)
                self.time_grant -= step

        if not self.applying and not self.made_redundant:
            self.check_reentry(params, rng)

        self.check_time_bounds()

    def update_strategy_self_thermostat(self, params, rng):

        """
        Update agent strategy based on self learning.

        1. agents who do not hold a grant increase their time allocation
        2. agents who hold a grant decrease their time allocation

        WCSS paper - THERMOSTAT model
        """

        step = _get_update_step(params, rng)

        if not self.grant_held:
            #self.time_grant *= (1.0 + step)
            self.time_grant += step
        else:
            #self.time_grant *= (1.0 - step)
            self.time_grant -= step

        if not self.applying:
            self.check_reentry(params, rng)

        self.check_time_bounds()


    # pylint: enable=too-many-instance-attributes

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
