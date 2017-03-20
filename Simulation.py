"""
--------------------------------------------------------------------------------
FILE: Simulation.py
DESC: Simulation functions
AUTH: thorsilver
VERS: 1.0 March 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""

from __future__ import print_function # ensure Python 3.x print form known
import pylab
from Population import Population
import Utils

class Simulation(object):

    """
    A single experimental run, together with data and derived statistics.
    """

    def __init__(self, params):

        self.params = params
        self.init_stats()
        self.sim = Population(params)
        self.colours = Utils.create_colour_dict(len(self.sim.agents))
        self.plot_colours = Utils.create_plot_colour_dict()
        self.endgame = int(self.params['iterations'] * 0.9)

    def init_stats(self):

        """
        Initialise data structures for storing data.
        """
        self.all_stats = []     # for endgame iterations
        self.final_stats = []   # for final iteration
        self.alloc = []
        self.accept = []

        self.all_tg = []
        self.tg_g = []
        self.tg_fg = []
        self.tg_ng = []

        self.all_r = []
        self.r_g = []
        self.r_fg = []
        self.r_ng = []
        self.r_pdr = []
        self.r_former_pdr = []
        self.r_old_academic = []

        self.all_rq = []
        self.rq_g = []
        self.rq_ng = []
        self.rq_fail = []
        self.rq_na = []
        self.rq_pdr = []
        self.rq_former_pdr = []
        self.rq_old_academic = []

        self.corr_rq_tg = []
        self.corr_rq_apply = []
        self.corr_rq_held = []

        self.pdr_total = []
        self.academic_total = []
        self.exited_total = []
        self.accepted_grants = []
        self.dynamic_roi = []
        self.roi_sum = []
        self.roi_sum_pdr = []


    def run(self):

        """
        Run the experiment.
        """
        iterations = 0
        for indx in range(self.params['iterations']): # *** was xrange ***
            iterations += 1
            self.sim.int_academic_count()
            self.sim.produce_applications()
            self.sim.evaluate_applications()
            self.sim.produce_research()
            if self.params['use_postdocs'] == 1:
                self.sim.hire_postdocs(self.params)
                self.sim.update_strategies()
            if self.params['growing_pop'] == 1:
                self.sim.update_newbies()
            if self.params['use_postdocs'] == 1:
                self.sim.update_postdocs()
            if self.params['use_retirement']:
                self.sim.update_careers()
                self.save_base_stats(t)
            if self.params['write_output']:
                self.save_plot_stats()
                print('Total agents: {}'.format(len(self.sim.agents)))
                self.sim.current_postdocs()
                self.accepted_grants.append(self.sim.acceptance_rate()[0])
                #print(self.accepted_grants)

            # self.dynamic_roi.append(
            #   self.calc_roi_dynamic(self.accepted_grants[indx])[2])
            # print('Current ROI result: {}'.format(self.dynamic_roi[indx]))
            self.roi_sum.append(self.calc_roi_sum(self.accepted_grants[indx]))
            #print('Current ROI sum: {}'.format(self.roi_sum[indx]))
            self.roi_sum_pdr.append(
                self.calc_roi_no_pdr(self.accepted_grants[indx]))
            #print('Current ROI, no PDRs: {}'.format(self.roi_sum_pdr[indx]))
            self.sim.clear_all()


    def test_flat(self, fixed_time):

        """
	Run a simulation in which all agents use a fixed time allocation.
	That is, no learning occurs.
	NB: there will still be some variation from iteration to iteration
	due to noise in the calculation of grant quality.
        """

        for agent in self.sim.agents:
            agent.time_grant = fixed_time

        # *** was xrange ***
        for _ in range(self.endgame, self.params['iterations']):
            self.sim.produce_applications()
            self.sim.evaluate_applications()
            self.sim.produce_research()
            self.save_base_stats(t)
            if self.params['write_output']:
                self.save_plot_stats()
            if self.params['use_postdocs'] == 1:
                self.sim.update_postdocs()
                self.sim.hire_postdocs()
            if self.params['growing_pop'] == 1:
                self.sim.update_newbies()

            self.sim.clear_all()


    def save_base_stats(self, thresh):

        """
        Save numerical stats related to efficiency.
        """

        if thresh >= self.endgame:
            self.all_stats.append(self.sim.all_stats())
        if thresh >= self.params['iterations'] - 1:
            self.final_stats = self.sim.all_stats()


    def save_plot_stats(self):

        """
        Save additional stats for plotting output.
        """

        alloc_accept = self.sim.acceptance_rate()
        self.alloc.append(alloc_accept[0])
        self.accept.append(alloc_accept[1])

        self.all_tg.append(self.sim.all_tg())
        self.tg_g.append(self.sim.all_tg_grant())
        self.tg_fg.append(self.sim.all_tg_fail())
        self.tg_ng.append(self.sim.all_tg_no_grant())

        self.all_r.append(self.sim.all_r())
        self.r_g.append(self.sim.all_r_grant())
        self.r_fg.append(self.sim.all_r_fail())
        self.r_ng.append(self.sim.all_r_no_grant())
        self.r_pdr.append(self.sim.all_r_pdr())
        self.r_former_pdr.append(self.sim.all_r_former_pdr())
        self.r_old_academic.append(self.sim.all_r_old_academic())

        self.all_rq.append(self.sim.all_rq())
        self.rq_g.append(self.sim.all_rq_grant())
        self.rq_ng.append(self.sim.all_rq_no_grant())
        self.rq_fail.append(self.sim.all_rq_fail())
        self.rq_na.append(self.sim.all_rq_no_apply())
        self.rq_pdr.append(self.sim.all_rq_pdr())
        self.rq_former_pdr.append(self.sim.all_rq_former_pdr())
        self.rq_old_academic.append(self.sim.all_rq_old_academic())

        self.pdr_total.append(self.sim.postdoc_count())
        self.academic_total.append(self.sim.academic_count())
        self.exited_total.append(self.sim.exited_count())

        self.corr_rq_tg.append(pylab.corrcoef(self.sim.all_rq(),
                                              self.sim.all_tg())[0][1])
        self.corr_rq_apply.append(pylab.corrcoef(self.sim.all_rq(),
                                                 self.sim.all_apply())[0][1])
        self.corr_rq_held.append(pylab.corrcoef(self.sim.all_rq(),
                                                self.sim.all_held())[0][1])


    def calc_mean_total_output(self):

        """
        *** Needs description ***
        """
        total_r_sum = 0
        count = 0
        for snap in self.all_stats:
            for a in snap:
                total_r_sum += a[5]
                count += 1
                #assert count == 0, "Div by zero, dipshit"
        if not count == 0:
            return total_r_sum / count
        else:
            return 0

    def calc_mean_time_grant(self):

        """
        *** Needs description ***
        """
        mean_tg_sum = 0
        count = 0
        for snap in self.all_stats:
            for a in snap:
                mean_tg_sum += a[3]
                count += 1
        return mean_tg_sum / count


    def calc_mean_corr_rq_held(self):

        """
        *** Needs description ***
        """
        corr_rq_held_sum = 0
        count = 0
        for snap in self.all_stats:
            rq_list = []
            held_list = []
            for a in snap:
                rq_list.append(a[1])
                held_list.append(a[4])
                corr_rq_held_sum += pylab.corrcoef(rq_list, held_list)[0][1]
                count += 1
        return corr_rq_held_sum / count


    def calc_roi(self, funding_units):

        """
        Calculate estimated return on investment obtained.
        """
        total_r = self.calc_mean_total_output()
        no_funding = self.sim.estimate_output(0.0, 0.0, 0.0, 'max')
        roi = (total_r - no_funding) / funding_units
        return total_r, no_funding, roi


    def calc_roi_dynamic(self, funding_units):

        """
        Calculated estimated ROI for dynamic population each timestep.
        """
        total_r = self.calc_mean_total_output()
        no_funding = self.sim.estimate_output(0.0, 0.0, 0.0, 'max')
        roi = ((total_r - no_funding)) / funding_units
        return total_r, no_funding, roi

    def calc_roi_sum(self, funding_units):

        """
        Calculate ROI - amount of extra total research bought per grant.
        """
        res_values = [a.research for agent in self.sim.agents
                      if not agent.made_redundant]
        total_res = sum(res_values)
        total_res_nof = self.sim.estimate_output_sum()
        roi = ((total_res - total_res_nof) - funding_units) / funding_units
        if len(res_values)>= 300 and roi >= -2.5:
            print('ROI jumps here: {} {} {}'.format(
                roi, total_res, total_res_nof))
            #quit()
        return roi


    def calc_roi_no_pdr(self, funding_units):

        """
        Calculate ROI - amount of research bought by employing PDRs.
        """

        res_values2 = [agent.research for agent in self.sim.agents
                       if not agent.made_redundant
                       and agent.postdoc_status == 0
                       and agent.former_postdoc == 0]
        total_res2 = sum(res_values2)
        total_res_nof2 = self.sim.estimate_output_sum()
        roi = ((total_res2 - total_res_nof2) - funding_units) / funding_units
        return roi

    def calc_redundancies(self):

        """
        Return total redundancies
        """
        self.redundancies_total = len(
            [agent for agent in self.sim.agents if agent.made_redundant])


    def write_output(self):

        """
        Write data plots.
        """
        # calculate means and sums of raw data
        t_all_tg = Utils.transpose(self.all_tg)
        mean_tg, sum_tg = Utils.get_mean_and_sum(self.all_tg)
        mean_tg_grant, sum_tg_grant = Utils.get_mean_and_sum(self.tg_g)
        mean_tg_fail, sum_tg_fail = Utils.get_mean_and_sum(self.tg_fg)
        mean_tg_no_grant, sum_tg_no_grant = Utils.get_mean_and_sum(self.tg_ng)

        t_all_r = Utils.transpose(self.all_r)
        self.mean_r, self.sum_r = Utils.get_mean_and_sum(self.all_r)
        mean_r_grant, sum_r_grant = Utils.get_mean_and_sum(self.r_g)
        mean_r_fail, sum_r_fail = Utils.get_mean_and_sum(self.r_fg)
        mean_r_no_grant, sum_r_no_grant = Utils.get_mean_and_sum(self.r_ng)
        self.mean_r_postdoc, sum_r_postdoc = Utils.get_mean_and_sum(self.r_pdr)
        self.mean_r_former_pdr, sum_r_former_pdr = Utils.get_mean_and_sum(
            self.r_former_pdr)
        self.mean_r_old_academic, sum_r_old_academic = Utils.get_mean_and_sum(
            self.r_old_academic)

        mean_rq, sum_rq = Utils.get_mean_and_sum(self.all_rq)
        mean_rq_grant, sum_rq_grant = Utils.get_mean_and_sum(self.rq_g)
        mean_rq_no_grant, sum_rq_no_grant = Utils.get_mean_and_sum(self.rq_ng)
        mean_rq_fail, sum_rq_fail = Utils.get_mean_and_sum(self.rq_fail)
        mean_rq_na, sum_rq_na = Utils.get_mean_and_sum(self.rq_na)
        mean_rq_pdr, sum_rq_pdr = Utils.get_mean_and_sum(self.rq_pdr)
        mean_rq_former_pdr, sum_rq_former_pdr = Utils.get_mean_and_sum(
            self.rq_former_pdr)
        mean_rq_old_academic, sum_rq_old_academic = Utils.get_mean_and_sum(
            self.rq_old_academic)

        self.redundancies_total = len(
            [agent for agent in self.sim.agents if agent.made_redundant])

        Utils.write_data_2d(self.sim.funding_body.successful_app_stats,
                            self.params['prefix']+"successful_app_stats.csv")

        # save final iteration stats
        Utils.write_data_2d(self.final_stats,
                            self.params['prefix']+"final_stats.csv")


        #### WRITE PLOTS

        # number of grants vs research quality
        Utils.write_plot([self.sim.all_rq()], [self.sim.all_grant_counts()],
                         self.params['prefix']+"lifetime_grants", "lifetime grants",
                         'research quality', 'grants awarded', self.plot_colours,
                         ylim = (0, int(self.params['iterations'])), marker='o',lw=0)

        # number of grants awarded
        #Utils.write_plot([range(len(self.alloc))], [self.alloc],
        #       self.params['prefix']+"grants_awarded", "grants awarded",
        #      'iterations', 'grants_awarded', self.plot_colours,
        #     ylim = (0, int(self.params['grant_proportion'] *
        #        self.params['pop_size']))

        # total postdocs

        Utils.write_plot([range(len(self.pdr_total))], [self.pdr_total],
                         self.params['prefix']+"total_pdrs", "Total Postdocs",
                         'Iterations', 'Total postdocs', self.plot_colours)

        # total postdocs vs total academics

        Utils.write_plot(
            [range(len(self.pdr_total)), range(len(
                self.academic_total)), range(len(self.exited_total))],
            [self.pdr_total, self.academic_total, self.exited_total],
            self.params['prefix']+"Academics vs Postdocs", '',
            'Iterations', 'Population', self.plot_colours,
            labels=('Postdocs', 'Academics', 'Leavers'))

        # acceptance rate (allocated / submitted)
        Utils.write_plot([range(len(self.accept))], [self.accept],
                         self.params['prefix']+"accept_rate", "acceptance rate",
                         'iterations', 'acceptance rate', self.plot_colours,
                         ylim = (0, 1))

        # # all time_grant values
        # Utils.write_plot([range(len(t_all_tg[0]))]*len(t_all_tg), t_all_tg,
        #         self.params['prefix']+"tg_all", "time_grant",
        #         'iterations', 'time_grant', self.colours,
        #         keys = self.sim.all_rq())

        # all research values
        # Utils.write_plot([range(len(t_all_r[0]))]*len(t_all_r), t_all_r,
        #         self.params['prefix']+"r_all", "research",
        #         'iterations', 'research', self.colours,
        #         keys = self.sim.all_rq())

        # mean time_grant
        Utils.write_plot(
            [range(len(mean_tg)), range(len(mean_tg_grant)),
             range(len(mean_tg_no_grant)), range(len(mean_tg_fail))],
            [mean_tg, mean_tg_grant, mean_tg_no_grant, mean_tg_fail],
            self.params['prefix']+"tg_mean", '', 'iterations',
            'mean time_grant', self.plot_colours,
            labels=('all', 'grant', 'no grant', 'fail'))


        #ROI EXPERIMENT 1
        # Utils.write_plot([range(len(self.dynamic_roi))], [self.dynamic_roi],
        #         self.params['prefix']+"ROI_test",
        #         "Return on Investment",
        #         'Iterations', 'ROI', self.plot_colours)

        #ROI EXPERIMENT 2
        Utils.write_plot(
            [range(len(self.roi_sum))], [self.roi_sum],
            self.params['prefix']+"ROI_test_sum", "Return on Investment",
            'Iterations', 'ROI', self.plot_colours)

        # mean research - no postdocs

        Utils.write_plot(
            [range(len(self.mean_r)), range(len(mean_r_grant)),
             range(len(mean_r_no_grant)), range(len(mean_r_fail))],
            [self.mean_r, mean_r_grant, mean_r_no_grant, mean_r_fail],
            self.params['prefix']+"r_mean", '', 'Iterations',
            'Mean Research Productivity', self.plot_colours,
            labels=('All', 'Grant', 'No Grant', 'Fail'))

        # mean research - postdocs
        Utils.write_plot(
            [range(len(self.mean_r)), range(len(mean_r_grant)),
             range(len(mean_r_no_grant)), range(len(mean_r_fail)),
             range(len(self.mean_r_postdoc))],
            [self.mean_r, mean_r_grant, mean_r_no_grant, mean_r_fail,
             self.mean_r_postdoc], self.params['prefix']+"r_mean_pdr", '',
            'Iterations', 'Mean Research Productivity', self.plot_colours,
            labels=('All', 'Grant', 'No Grant', 'Fail', 'PDRs'))

        # mean research - postdocs, variant set
        Utils.write_plot(
            [range(len(self.mean_r)), range(len(mean_r_grant)),
             range(len(mean_r_no_grant)), range(len(self.mean_r_former_pdr)),
             range(len(self.mean_r_postdoc))],
            [self.mean_r, mean_r_grant, mean_r_no_grant,
             self.mean_r_former_pdr, self.mean_r_postdoc],
            self.params['prefix']+"r_mean_pdr2", '', 'Iterations',
            'Mean Research Productivity', self.plot_colours,
            labels=('All', 'Grant', 'No Grant', 'Promoted PDRs', 'PDRs'))

        # mean research - academics vs promoted postdocs
        Utils.write_plot(
            [range(len(self.mean_r)), range(len(self.mean_r_postdoc)),
             range(len(self.mean_r_former_pdr)),
             range(len(self.mean_r_old_academic))],
            [self.mean_r, self.mean_r_postdoc, self.mean_r_former_pdr,
             self.mean_r_old_academic],
            self.params['prefix']+"r_mean_promoted_pdr", '', 'Iterations',
            'Mean Research Productivity', self.plot_colours,
            labels=('All', 'PDRs', 'Promoted PDRs', 'Established Academics'))

        # sum time_grant
        Utils.write_plot(
            [range(len(sum_tg)), range(len(sum_tg_grant)),
             range(len(sum_tg_no_grant)), range(len(sum_tg_fail))],
            [sum_tg, sum_tg_grant, sum_tg_no_grant, sum_tg_fail],
            self.params['prefix']+"tg_sum", '',
            'iterations', 'sum time_grant', self.plot_colours,
            labels=('all', 'grant', 'no grant', 'fail'))

        # sum research
        Utils.write_plot(
            [range(len(self.sum_r)), range(len(sum_r_grant)),
             range(len(sum_r_no_grant)), range(len(sum_r_fail))],
            [self.sum_r, sum_r_grant, sum_r_no_grant, sum_r_fail],
            self.params['prefix']+"r_sum", '',
            'iterations', 'sum research', self.plot_colours,
            labels=('all', 'grant', 'no grant', 'fail'))

        # sum research - postdocs
        Utils.write_plot(
            [range(len(self.sum_r)), range(len(sum_r_grant)),
             range(len(sum_r_no_grant)), range(len(sum_r_fail)),
             range(len(sum_r_postdoc))],
            [self.sum_r, sum_r_grant, sum_r_no_grant, sum_r_fail,
             sum_r_postdoc], self.params['prefix']+"r_sum_pdr", '',
            'Iterations', 'Sum Research', self.plot_colours,
            labels=('All', 'Grant', 'No Grant', 'Fail', 'Postdoc'))

        # sum research - promoted postdocs vs established academics
        Utils.write_plot(
            [range(len(self.sum_r)), range(len(sum_r_postdoc)),
             range(len(sum_r_former_pdr)), range(len(sum_r_old_academic))],
            [self.sum_r, sum_r_postdoc, sum_r_former_pdr, sum_r_old_academic],
            self.params['prefix']+"r_sum_promoted_pdr", '',
            'Iterations', 'Sum Research', self.plot_colours,
            labels=('All', 'PDRs', 'Promoted PDRs', 'Established Academics'))

        # mean research_quality of successful applicants - no postdocs
        Utils.write_plot(
            [range(len(mean_rq)), range(len(mean_rq_grant)),
             range(len(mean_rq_no_grant)), range(len(mean_rq_fail)),
             range(len(mean_rq_na))],
            [mean_rq, mean_rq_grant, mean_rq_no_grant, mean_rq_fail, mean_rq_na],
            self.params['prefix']+"x_rq_spread", '',
            'Iterations', 'Mean Research Quality', self.plot_colours,
            labels=('All', 'Grant', 'No Grant', 'Fail', 'Not Apply'))

        # mean research_quality of successful applicants - postdocs
        Utils.write_plot(
            [range(len(mean_rq)), range(len(mean_rq_grant)),
             range(len(mean_rq_no_grant)),
             range(len(mean_rq_fail)), range(len(mean_rq_pdr))],
            [mean_rq, mean_rq_grant, mean_rq_no_grant, mean_rq_fail,
             mean_rq_pdr], self.params['prefix']+"x_rq_pdr", '', 'Iterations',
            'Mean Research Quality', self.plot_colours,
            labels=('All', 'Grant', 'No Grant', 'Fail', 'Postdoc'))

        # correlations: rq v tg
        Utils.write_plot(
            [range(len(self.corr_rq_tg)), range(len(self.corr_rq_apply)),
             range(len(self.corr_rq_held))],
            [self.corr_rq_tg, self.corr_rq_apply, self.corr_rq_held],
            self.params['prefix']+"corr_rq", '',
            'Iterations', 'Correlation Coef.', self.plot_colours,
            labels = ('(RQ, TG)', '(RQ, Apply)', '(RQ, Held)'),
            ylim = (-0.75,0.75))
