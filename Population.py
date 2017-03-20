"""
--------------------------------------------------------------------------------
FILE: Population.py
DESC: ***
AUTH: thorsilver
VERS: 1.0 March 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""

# *** Unused imports ***
# import collections
# import operator

from __future__ import print_function # ensure Python 3.x print form known
from random import Random
from random import randint
import random
import math
from FundingAgency import FundingAgency
from Academic import Academic
from Application import Application


class Population(object):

    """
    A population of academic agents.

    This class handles population level actions such as producing grants and
    research, and allocating funding.  It also provides an interface to obtain
    data/statistics on agents.
    """

    def __init__(self, params):
        self.params = params
        self.postdocs_generated = 0
        self.timestamp = 0
        self.leavers = 0
        self.hired = 0
        self.academics = 0
        self.new_id = 0         # *** added by IW - needs checking ***
        self.postdoc_id = 0     # *** added by IW - needs checking ***

        if self.params['random_seed']:
            self.rng = Random()
        else:
            self.rng = Random(self.params['seed'])

        # initialise agents
        self.funding_body = FundingAgency(params)
        self.agents = []
        for i in range(self.params['pop_size']): # *** was xrange ***
            self.agents.append(Academic(i, params, self.rng))

        # calculate derived parameters
        params['grant_count'] = int(
            len(self.agents) * self.params['grant_proportion'])


    def estimate_output(self, bonus, prop, time=0.0, typ='rnd'):

        """
        Estimate the total research output of a population given:
        - typ='max' : best possible allocation of grants
        - typ='rnd' : random allocation of grants
			(averaged over several attempts)
        - typ='min' : worst possible allocation of grants
        and that each individual spends a fixed and equal amount of
        time on their applications.
        """

        attempts = 1

        rq_agents = [(agent.research_quality, agent)
                     for agent in self.agents if not agent.made_redundant]
        if typ == 'max':
            rq_agents.sort(reverse=True)
        elif typ == 'min':
            rq_agents.sort()
        elif typ == 'rnd':
            attempts = 10

        research_sum = 0.0
        grant_number = int(len(rq_agents) * prop)

        for _ in range(attempts):
            for agent in rq_agents[:grant_number]:
                research_sum += agent[1].calc_research(
                    time, True, bonus, agent[0])
            for agent in rq_agents[grant_number:]:
                research_sum += agent[1].calc_research(
                    time, False, bonus, agent[0])

        return research_sum / attempts


    def estimate_output_sum(self):

        """
        Estimate total research across system without any grant funding.
        """
        research_sum = 0.0
        rq_agents = [(agent.research_quality, agent)
                     for agent in self.agents if not agent.made_redundant]
        rq_agents.sort(reverse=True)
        for agent in rq_agents:
            research_sum += agent[1].calc_research(0.0, True, 0.0, agent[0])
        return research_sum


    ## SIM ACTIONS ############################################################

    def produce_applications(self):

        """
        Produce applications by each agent (who is applying).
        """
        [self.funding_body.add_application(
            Application(agent, self.params, self.rng), self.rng)
         for agent in self.agents if agent.applying]


    def evaluate_applications(self):

        """
        Evaluate the submitted applications and allocate grants. Generate
        postdocs in the postdoc case, or new academics in the simple
        growing population case

        """
#        if self.params['limited_funding'] == True and iterations == 1:
#            self.funding_body.init_grants()
        postdoc_noise = self.rng.randint(1, 4)
        grant_total = 0
        self.funding_body.rank_applications()
        if not self.params['limited_funding']:
            successful = self.funding_body.get_grant_recipients(
                self.params, self.academics)
        else:
            successful = self.funding_body.get_recipients_limited(
                self.params, self.academics)
        print("Total academics now: {}".format(self.academics))
        # print(successful)
        for agent in successful:
            self.agents[agent].grant_held = True
            self.agents[agent].grant_count += 1
            if self.params['use_postdocs'] == 1 or \
               self.params['growing_pop'] == 1:
                grant_total += 1
        if self.params['use_postdocs'] == 1:
            self.generate_postdocs(grant_total//10 + postdoc_noise)
        if self.params['growing_pop'] == 1:
            self.add_academics(grant_total//20 + postdoc_noise)

        print('Total grants disbursed: {}'.format(grant_total))




##    def evaluate_applications_pdr(self):
##
##        """
##        Evalute the submitted applications and allocate grants.
##        """
##
##        postdoc_noise = self.rng.randint(1, 4)
##        grant_total = 0
##        self.funding_body.rank_applications()
##        successful = self.funding_body.get_grant_recipients_pdr(
##             self.params, len(self.agents) - self.postdocs_generated)
##        for agent in successful:
##            self.agents[a].grant_held = True
##            self.agents[a].grant_count += 1
##            if self.params['use_postdocs'] == 1:
##                grant_total += 1
##        if self.params['use_postdocs'] == 1:
##            self.generate_postdocs(grant_total//6 + postdoc_noise)
##        print("Total cash moneyz: {}".format(grant_total))


    def produce_research(self):

        """
        Produce research by each agent.  Return total research.
        """

        return sum([agent.produce_research(self.params)
                    for agent in self.agents if not agent.made_redundant])


    def update_strategies(self):

        """
        Update agent strategies.
        """
        for agent in self.agents:
            if agent.postdoc_status == 0:
                if self.params['learning_type'] == 'thermostat':
                    agent.update_strategy_self_thermostat(self.params, self.rng)
                elif self.params['learning_type'] == 'memory':
                    agent.update_strategy_self_memory(self.params, self.rng)
                else:
                    System.exit("Unknown learning type")


    def clear_all(self):

        """
        Clear any grants currently held by agents.
        """
        for agent in self.agents:
            agent.grant_held = False
        self.funding_body.clear_applications()

    def update_postdocs(self):

        """
        Update postdoc contracts and status
        """
        for agent in self.agents:
            if agent.postdoc_status == 1:
                if agent.contract_length >= 1:
                    agent.contract_length -= 1
                if agent.newb >= 1:
                    agent.newb -= 1
                # print('Contract Length: {}'.format(agent.contract_length))

    def update_careers(self):

        """
        Decrement career_length by one each iteration
        """
        for agent in self.agents:
            if agent.career_length >= 1:
                agent.career_length -= 1
            elif agent.career_length <= 1 and not agent.retired:
                agent.retire()

    def update_newbies(self):
        """
        Update newbie academics in growing population case
        """
        for agent in self.agents:
            if agent.newb >= 1:
                agent.newb -= 1



    def hire_postdocs(self, params):

        """
        Base: Hire 15% of postdocs (lucky!)
        RQ counts: postdocs with higher RQ get hired
        """

        leavers = [agent for agent in self.agents
                   if agent.postdoc_status == 1 and agent.contract_length <= 0
                   and not agent.made_redundant]
        promotions_count = 0
        redundancies_count = 0
        if len(leavers) > 0:
            ranked_leavers = sorted(
                leavers, key=lambda x: x.research_quality, reverse=True)
            candidate_num = int(math.ceil(
                params['postdoc_chance'] * len(leavers)))
            #print('Research Quality counts: {}'.format(
            #  self.params['pdr_rq_counts']))
            #print('Postdoc Chance: {}'.format(params['postdoc_chance']))
            #print('Candidates: {}'.format(candidate_num))
            #print('Ranked Leavers: {}'.format(len(ranked_leavers)))
            if self.params['pdr_rq_counts'] == 0:
                random.shuffle(leavers)
                for agent in ranked_leavers[:candidate_num]:
                    agent.former_postdoc = 1
                    agent.tenured = True
                    agent.postdoc_status = 0
                    agent.set_random_time_grant(
                        self.params['init_time'], self.rng)
                    agent.applying = True
                    promotions_count += 1
                for agent in ranked_leavers[candidate_num:]:
                    agent.made_redundant = True
                    agent.exited_postdoc = 1
                    redundancies_count += 1
            elif self.params['pdr_rq_counts'] == 1:
                for agent in ranked_leavers[:candidate_num]:
                    agent.former_postdoc = 1
                    agent.tenured = True
                    agent.postdoc_status = 0
                    agent.set_random_time_grant(self.params['init_time'], self.rng)
                    agent.applying = True
                    promotions_count += 1
                for agent in ranked_leavers[candidate_num:]:
                    agent.made_redundant = True
                    agent.exited_postdoc = 1
                    redundancies_count += 1
            #print('Total agents sacked this round: {}'.format(redundancies_count))
            #print('Total agents promoted this round: {}'.format(promotions_count))
            if self.params['pdr_rq_counts'] == 0:
                self.agents = list(set(self.agents+leavers))
            else:
                self.agents = list(set(self.agents + ranked_leavers))


    def add_academics(self, num_academics):

        """
        Add new academics for agent simple growing population case
        Each timestep a selection of new academics come in based on the number of grants issued
        """
        for _ in range(0, num_academics):
            self.new_id = (len(self.agents))
            print("Length of agent list: {}".format(self.new_id))
            self.agents.append(Academic(self.new_id, self.params, self.rng))
            self.agents[self.new_id].research_quality = self.rng.random()
            self.agents[self.new_id].applying = True
            self.agents[self.new_id].grant_held = False
            self.agents[self.new_id].tenured = True
            self.agents[self.new_id].grant_count = 0
            self.agents[self.new_id].set_random_time_grant(
                self.params['init_time'], self.rng)
            self.agents[self.new_id].research = 0.0
            self.agents[self.new_id].research_sum = 0.0
            self.agents[self.new_id].memory = []
            self.agents[self.new_id].former_postdoc = 0
            self.agents[self.new_id].exited_postdoc = 0
            self.agents[self.new_id].num_postdocs = 0
            self.agents[self.new_id].postdoc_status = 0
            self.agents[self.new_id].contract_length = 100
            self.agents[self.new_id].time_grant = 0.0
            self.agents[self.new_id].newb = 2


    def generate_postdocs(self, num_postdocs):

        """
        Generate postdocs as requested by grant allocation functions
        Each new postdoc gets assigned appropriate attributes
        """
        for _ in range(0, num_postdocs):
            self.postdoc_id = (len(self.agents))
            #print("New postdoc: {}".format(self.postdoc_id))
            self.agents.append(Academic(
                self.postdoc_id, self.params, self.rng))
            self.agents[self.postdoc_id].research_quality = self.rng.random()
            self.agents[self.postdoc_id].applying = False
            self.agents[self.postdoc_id].grant_held = False
            self.agents[self.postdoc_id].tenured = False
            self.agents[self.postdoc_id].grant_count = 0
            self.agents[self.postdoc_id].set_random_time_grant(
                self.params['init_time'], self.rng)
            self.agents[self.postdoc_id].research = 0.0
            self.agents[self.postdoc_id].research_sum = 0.0
            self.agents[self.postdoc_id].memory = []
            self.agents[self.postdoc_id].former_postdoc = 0
            self.agents[self.postdoc_id].exited_postdoc = 0
            self.agents[self.postdoc_id].num_postdocs = 1
            self.agents[self.postdoc_id].postdoc_status = 1
            self.agents[self.postdoc_id].contract_length = randint(4, 10)
            self.agents[self.postdoc_id].time_grant = 0.0
            self.agents[self.postdoc_id].newb = 2

##            self.agents[a].research_quality = self.rng.random()
##            self.agents[a].applying = False
##            self.agents[a].grant_held = False
##            self.agents[a].tenured = False
##            self.agents[a].grant_count = 0
##            self.agents[a].set_random_time_grant(self.params['init_time'], self.rng)
##            self.agents[a].research = 0.0
##            self.agents[a].research_sum = 0.0
##            self.agents[a].memory = []
##            self.agents[a].former_postdoc = 0
##            self.agents[a].exited_postdoc = 0
##            self.agents[a].num_postdocs = 1
##            self.agents[a].postdoc_status = 1
##            self.agents[a].contract_length = randint(4,10)
##            self.agents[a].time_grant = 0.0
##            self.agents[a].newb = 2
            self.postdocs_generated += 1


    ## DATA ACCESS ############################################################

    def all_stats(self):

        """
        Return a table of (for the current iteration):
        [ID, rq, app, tg, g, r, PDR status]
        """
        return [
            (agent.id, agent.research_quality, agent.applying,
             agent.time_grant, agent.grant_held, agent.research,
             agent.postdoc_status) for agent in self.agents]

    def acceptance_rate(self):

        """
        Return tuple containing # grants allocated and acceptance rate:
        (# grants allocated, (# grants allocated) / (# grants submitted)).
        """
        submitted = 0
        allocated = 0
        for agent in self.agents:
            if agent.applying:
                submitted += 1
            if agent.grant_held:
                allocated += 1
        if submitted > 0:
            return allocated, float(allocated) / submitted
        else:
            return allocated, 0.0


    def all_rq(self):

        """
        Return a list of all research quality values.
        """
        return [agent.research_quality for agent in self.agents]


    def all_r(self):

        """
        Return a list of all research output values for current year.
        """
        return [agent.research for agent in self.agents]


    def all_tg(self):

        """
        Return a list of all tg values.
        """
        return [agent.time_grant for agent in self.agents]


    def all_apply(self):

        """
        Return a list of all applying values.
        """
        return [agent.applying for agent in self.agents]


    def all_held(self):

        """
        Return a list of all grant_held values.
        """
        return [agent.grant_held for agent in self.agents]


    def all_r_grant(self):

        """
        Return a list of tg values of agents holding grants.
        """
        return [agent.research for agent in self.agents if agent.grant_held]


    def all_r_fail(self):

        """
        Return a list of r values of agents who apply but fail.
        """
        return [agent.research for agent in self.agents if
                (agent.applying and not agent.grant_held)]


    def all_r_no_grant(self):

        """
        Return a list of r values of agents not holding grants.
        """
        return [agent.research for agent in self.agents if not agent.grant_held]


    def all_r_pdr(self):

        """
        Return a list of r values of agents who are postdocs.
        """
        return [agent.research for agent in self.agents
                if agent.postdoc_status == 1 and not agent.made_redundant]

    def all_r_former_pdr(self):

        """
        Return a list of r values of agents who were promoted.
        """
        return [agent.research for agent in self.agents
                if agent.former_postdoc == 1 and not agent.made_redundant]


    def all_r_old_academic(self):

        """
        Return a list of r values of agents who are established academics.
        """
        #return [agent.research for agent in self.agents if agent.postdoc_status == 0
        # and agent.former_postdoc == 0 and not agent.made_redundant]
        return [agent.research for agent in self.agents
                if agent.postdoc_status == 0 and agent.former_postdoc == 0
                and agent.career_length <= 20]

    def all_tg_grant(self):

        """
        Return a list of tg values of agents holding grants.
        """
        return [agent.time_grant for agent in self.agents
                if agent.grant_held]


    def all_tg_fail(self):

        """
        Returns a list of tg values of agents who apply but fail.
        """
        return [agent.time_grant for agent in self.agents if
                (agent.applying and not agent.grant_held)]


    def all_tg_no_grant(self):

        """
        Return a list of tg values of agents not holding grants.
        """
        return [agent.time_grant for agent in self.agents
                if not agent.grant_held]


    # *** already defined earlier in this code ***
    # def all_rq(self):

    #     """
    #     Return a list of rq values of all agents.
    #     """
    #     return [agent.research_quality for agent in self.agents]


    def all_rq_grant(self):

        """
        Return a list of rq values of agents holding grants.
        """
        return [agent.research_quality for agent in self.agents
                if agent.grant_held]


    def all_rq_no_grant(self):

        """
        Return a list of rq values of agents not holding grants.
        """
        return [agent.research_quality for agent in self.agents
                if not agent.grant_held]


    def all_rq_fail(self):

        """
        Returns a list of rq values of agents who apply but fail.
        """
        return [agent.research_quality for agent in self.agents if
                (agent.applying and not agent.grant_held)]


    def all_rq_no_apply(self):

        """
        Returns a list of rq values of agents who don't apply.
        """
        return [agent.research_quality for agent in self.agents
                if not agent.applying]


    def all_rq_pdr(self):

        """
        Returns a list of rq values of agents who are postdocs.
        """
        return [agent.research_quality for agent in self.agents
                if agent.postdoc_status == 1 and not agent.made_redundant]

    def all_rq_former_pdr(self):

        """
        Returns a list of RQ values of agents who have been promoted.
        """
        return [agent.research_quality for agent in self.agents
                if agent.former_postdoc == 1 and not agent.made_redundant]

    def all_rq_old_academic(self):

        """
        Returns a list of RQ values of agents who are established academics.
        """
        return [agent.research_quality for agent in self.agents
                if agent.postdoc_status == 0 and agent.former_postdoc == 0
                and not agent.made_redundant]


    def all_grant_counts(self):

        """
        Returns a list of lifetime grant counts for each agent.
        """
        return [agent.grant_count for agent in self.agents]


    def postdoc_count(self):

        """
        Returns a total population of postdocs ever generated.
        """
        return [self.postdocs_generated]


    def academic_count(self):

        """
        Returns a total number of academics.
        """
        academics = 0
        for agent in self.agents:
            if agent.tenured:
                academics += 1
        return [academics]


    def int_academic_count(self):

        """
        Adds academic count to agent list.
        """

        self.academics = 0
        for agent in self.agents:
            if agent.tenured:
                self.academics += 1
        return [self.academics]

##        tenured = 0
##        new_faculty = 0
##        leavers = 0
##        for agent in self.agents:
##            if agent.former_postdoc == 1:
##                new_faculty += 1
##        for agent in self.agents:
##            if agent.exited_postdoc == 1:
##                leavers += 1
##        for agent in self.agents:
##            if agent.applying:
##                tenured += 1
##        return [tenured - self.postdocs_generated + new_faculty - leavers]

    def exited_count(self):

        """
        Returns number of exited postdocs.
        """
        self.leavers = 0
        for agent in self.agents:
            if agent.made_redundant:
                self.leavers += 1
        print('Leavers: {}"  Total Postdocs: {}'.format(
            self.leavers, self.postdocs_generated))
        return [self.leavers]


    def current_postdocs(self):

        """
        Returns number of active postdocs in the current semester.
        """
        active_pdrs = 0
        pdrs = 0
        for agent in self.agents:
            if agent.postdoc_status == 1 and agent.former_postdoc == 0 \
            and not agent.made_redundant:
                active_pdrs += 1
        for agent in self.agents:
            if agent.postdoc_status == 1:
                pdrs += 1
        print('Active Postdocs:  {}  Other postdocs:  {}'.format(
            active_pdrs, pdrs))
        # return [active_pdrs]


    def print_all(self):

        """
        *** do something ***
        """
        for agent in self.agents:
            agent.print_all()

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
