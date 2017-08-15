"""
------------------------------------------------------------------------
FILE: Dispatch.py
DESC: Dispatch simulations
AUTH: thorsilver
UPDT: digmusic
VERS: 1.1 July 2017
REQS: Python 3.x (version 3.6 used)
------------------------------------------------------------------------
"""

from __future__ import print_function # ensure Python 3.x print form known
import copy                           # for dictionary copying
import json                           # for JSON init file handling
import math
import os                             # operating system support 
import pylab
import random                         # random number generation 

# Installed via Macports
import matplotlib.pyplot as plt
import pandas as pd
# MAJS files
from Simulation import Simulation
import Utils

# Location of preferences file - will be different for a Windows-based
# implementation (need to sort this out). Contents of this file
# leads to where test run data files can be found.
PREFERENCES_FILE = os.path.expanduser("~/.majs")

# Filename constants with omitted extensions
INCRUN_MEAN = "incRunSet_meanR0"
INCRUN_ROI = "incRunSet_ROI"
INCRUN_SACK = "incRunSet_sacked"
INCRUN_TOT = "incRunSet_totalRO"
MEMRUN_MEAN = "memRunSet_meanRO"
MEMRUN_ROI = "memRunSet_ROI"
MEMRUN_SACK = "memRunSet_sacked"
MEMRUN_TOT = "memRunSet_totalRO"
PRORUN_MEAN = "promoRunSet_meanRO"
PRORUN_ROI = "promoRunSet_ROI"
PRORUN_SACK = "promoRunSet_sacked"
PRORUN_TOT = "promoRunSet_totalRO"

FUNDINC = 'Funding_Increase'
FUNDLIM = 'Funding_Limit'
GROWPOP = 'Growing_Pop'
MEANROALL = 'Mean_RO_all'
PROCHAN = 'Promo_Chance'
ROINOPDRS = 'ROI_no_PDRs'
TOTRO = 'Total_RO'
TOTSACK = 'Total_Sacked'
USEPDRS = 'Use_PDRs'

# ============================================================
# Private helper functions
# ============================================================

def _headers(dfile, stress, mean, sack, fund):
    """
    Write headers.
    If any of the arguments after the data file spec are
    true, the appropriate header part is written
    """
    dfile.write("Run Number,")
    if stress:
        dfile.write("Job_Stress,")
    dfile.write("Use_PDRs,Growing_Pop,RQ_counts,Mentored_PDRs,ROI," + \
                "ROI no PDRs,Total RO,")
    if mean:
        dfile.write("Mean RO all,")
    dfile.write("Mean RO Old Farts,Mean RO PDR,Mean RO FPDR,Mean TG," + \
                "Learning Type,Promo Chance")
    if sack:
        dfile.write(",Total Sacked")
    if fund:
        dfile.write(",Funding Limit,Funding Increase")
    dfile.write("\n")


def _init_parms(parms_file):
    """
    Initialise parameters at start of run. Normally values are
    taken from the default JSON file: "json/_default.json".
    See __init__ for the call.
    """

    with open(parms_file) as json_data:
        parms = json.load(json_data)

    return parms


def _init_prefs(prefs_file):
    """
    Initialise preferences from the .majs file in users' home directory
    """

    with open(prefs_file) as json_data:
        prefs = json.load(json_data)

    return prefs


def _plot_fig(depfig, item, ptitle, errval, fcolor, xlabel, ylabel,
              ticklbl, rootdir, fname, clflag=True):
    """
    Plot a figure
    """
    if clflag:
        plt.clf()
    fig1 = item.plot(
        kind='bar', title=ptitle, yerr=errval, facecolor=fcolor,
        alpha=0.5, rot=0)
    fig1.set_xlabel(xlabel, fontsize=12)
    fig1.set_ylabel(ylabel, fontsize=12)
    fig1.set_xticklabels(ticklbl)
    if depfig:
        fig2 = depfig.get_figure()
    else:
        fig2 = fig1.get_figure()
    fig2.savefig(rootdir + fname + '.pdf')
    fig2.savefig(rootdir + fname + '.png')
    return fig1


def _tabs(strs, nline=True):
    """
    Construct a string of the strs argument, separated by tabs. The
    default is to end the string with a newline - passing a second arg
    of False will omit the newline
    """

    last_item = strs.pop()
    ret = ""
    for item in strs:
        ret += str(item) + "\t" # str ensures a string even if already so
    ret += last_item
    if nline:
        ret += "\n"

    return ret

# ------------------------------------------------------------
# Dispatch class
# ------------------------------------------------------------

class Dispatch(object):
    """
    Groups together the simulation methods and preferences/parameters
    they require
    """

    def __init__(self):
        """
        Initialise parameters and preferences. Note that the preferences
        file defines the default file for parameters
        """

        self.prefs = _init_prefs(os.path.expanduser(PREFERENCES_FILE))
        # Note that the values in self.parms should not be overwritten
        # as they are the default parameters.
        # Instead use myparms = {**self.parms, **new_parms} to create a
        # local copy of the parms overwritten by those values in
        # new_parms, where parms and new_parms are dictionaries
        self.parms = _init_parms(self.prefs['default_parms'])

    # --------------------------------------------------
    # Private Dispatch methods - do not call directly
    # --------------------------------------------------

    def _load_data_file(self, uniqname):
        """
        Load a JSON data file. Uniqname e.g. "sim1" translates into
        the fully qualified name "json/sim1.json". The prefs dictionary
        contains the directory prefix and file extension values
        """

        fname = self.prefs['json_subdir'] + uniqname + \
                self.prefs['json_file_ext']
        print(fname)
        with open(fname) as json_data:
            return json.load(json_data)


    def _run_sim(self, dfile, parms, exp_prefix, sim=True):
        """
        FUNC: Run multiple initialisations of single set of experimental
        conditions.
        ARGS: dfile      - data file
              parms      - parameters to control the runs
              exp_prefix - prefix for experiments
              sim        - simulation flag. If sim=False, output is
                           estimated, rather than simulated.  Currently,
                           this is just used for evaluating the "fixed
                           time investment" heuristic.
        RETS:
        """
        sum_roi_values = []
        roi_values = []
        total_r_values = []
        mean_tg_values = []
        corr_rq_held_values = []
        redundancies_final = 0
        lparms = copy.deepcopy(parms) # local copy of parameters arg

        # for each initialisation
        for run in range(lparms['runs']):
            lparms['seed'] = random.randint(0, 99999999)
            lparms['prefix'] = exp_prefix + "%d/" % lparms['seed']
            print("   run #%d" % run)

            # run experiment; either simulation or estimation
            exp = Simulation(lparms)
            if sim:
                threshold = 10 # *** threshold needs looking at ***
                exp.run(threshold)
            else:
                exp.test_flat(lparms['fixed_time'])

            roi_values.append(exp.calc_roi(1)[2])
            total_r_values.append(exp.calc_mean_total_output())
            mean_tg_values.append(exp.calc_mean_time_grant())
            corr_rq_held_values.append(exp.calc_mean_corr_rq_held())
            sum_roi_values.append(exp.calc_roi_sum)
            #roi_SEs.append(numpy.std(sum_roi_values))
            exp.calc_redundancies()
            redundancies_final = exp.redundancies_total

            if lparms['write_output']:
                Utils.create_dir(lparms['prefix'])
                exp.write_output()
                #redundancies_final = exp.redundancies_total
                #redundancies_final.append(exp.redundancies_total)
                dfile.write(
                    str(run + 1) + "," + str(lparms['use_postdocs']) + \
                    "," + \
                    str(lparms['growing_pop']) + "," + \
                    str(lparms['pdr_rq_counts']) + "," + \
                    str(lparms['mentored_pdrs']) + "," + \
                    str(exp.roi_sum[lparms['iterations'] - 1]) + "," + \
                    str(exp.roi_sum_pdr[lparms['iterations'] - 1]) + \
                    "," + \
                    str(total_r_values[run]) + "," + \
                    str(exp.mean_r[lparms['iterations'] - 1]) + "," + \
                    str(exp.mean_r_old_academic[lparms['iterations'] - 1]) + \
                    "," + \
                    str(exp.mean_r_postdoc[lparms['iterations'] - 1]) + "," + \
                    str(exp.mean_r_former_pdr[lparms['iterations'] - 1]) + \
                    "," + \
                    str(mean_tg_values[run]) + "," + \
                    str(lparms['learning_type']) + "," + \
                    str(lparms['postdoc_chance']) + "," + \
                    str(exp.redundancies_total) + "," + \
                    str(lparms['limited_funding']) + "," + \
                    str(lparms['yearly_increase']) + "\n")

        # calculate and write summary statistics
        summ_data = [
            roi_values, total_r_values, mean_tg_values,
            corr_rq_held_values]

        summ_stats = [
            pylab.mean(roi_values), pylab.std(roi_values),
            pylab.mean(total_r_values), pylab.std(total_r_values),
            pylab.mean(mean_tg_values), pylab.std(mean_tg_values),
            pylab.mean(corr_rq_held_values),
            pylab.std(corr_rq_held_values)]

        if lparms['write_output']:
            Utils.write_data_2d(
                summ_data, exp_prefix + self.prefs['summ_data'] + \
                '.csv')
            Utils.write_data(
                summ_stats, exp_prefix + self.prefs['summ_stats'] + \
                '.csv')

        return pylab.mean(roi_values), pylab.std(roi_values), \
            total_r_values[-1], redundancies_final


    # ------------------------------------------------------------
    # Alife sim runs. Use func to indicate which sim to run
    # ------------------------------------------------------------

    def alife_run(self, func, seed, fund):
        """
        Run an Alife XV-based simulation.
        Expects the function to be called using parms, parms-prefix,
        and seed. Note: no need to copy parms as they are used as
        read-only
        """
        dfile = open(self.prefs['roi_file'], 'w')
        _headers(dfile, False, True, True, fund)
        rdir = self.parms['prefix']           # root directory
        func(self.parms, rdir, seed)
        dfile.close()
        return rdir

    # ------------------------------------------------------------
    # Alife sim runs
    # ------------------------------------------------------------

    # *** The code extracts parms that do not appear to change after
    #     initial value is set and places them outside of the loop.
    #     If this assumption is incorrect then the sims functions
    #     need to be looked at
    # ***

    def run_alife_sims(self, dfile, base_prefix, seed):
        """
        Run sims for Alife XV paper.
        """
        aparms = self._load_data_file('alife_sims')
        lparms = {**self.parms, **aparms}
        pdocs = 0
        indx = 0
        for simul in ["growingPop", "noRQnoM", "RQnoM", "RQM", "noRQM"]:
            lparms['growing_pop'] = 1 - pdocs
            lparms['mentored_pdrs'] = indx / 3
            lparms['pdr_rq_counts'] = (indx / 2) - (indx / 4)
            lparms['prefix'] = base_prefix + "memory_" + simul
            lparms['use_postdocs'] = pdocs
            print(simul)
            self._run_sim(dfile, lparms, lparms['prefix'])
            indx += 1
            if indx == 1:
                pdocs = 1 - pdocs

        return lparms

    # ------------------------------------------------------------
    # Alife sweeps 1
    # ------------------------------------------------------------

    def run_alife_sweep1(self, dfile, base_prefix):
        """
        Run sweep of promotion chance values.
        """
        aparms = self._load_data_file('alife_sweep1')
        lparms = {**self.parms, **aparms}
        for simul in lparms['loops']:
            pcent = str(int(100 * simul))
            lparms['postdoc_chance'] = simul
            lparms['prefix'] = base_prefix + "memory_chance" + pcent + "_"
            print("Memory, postdocs promo " + pcent + "percent")
            self._run_sim(dfile, lparms, lparms['prefix'])

        return lparms

    # ------------------------------------------------------------
    # Alife sweeps 2
    # ------------------------------------------------------------

    def run_alife_sweep2(self, dfile, base_prefix):
        """
        Run sweep of funding limits.
        """
        aparms = self._load_data_file('alife_sweep2')
        lparms1 = {**self.parms, **aparms['sim1']}
        print('Funding Limited')
        self._run_sim(dfile, lparms1, base_prefix)

        # run funding unlimited sim
        lparms2 = {**lparms1, **aparms['sim2']}
        print('Funding Unlimited')
        self._run_sim(dfile, lparms2, base_prefix)
        return lparms2

    # ------------------------------------------------------------
    # Alife sweeps 3
    # ------------------------------------------------------------

    def run_alife_sweep3(self, dfile, base_prefix):
        """
        Run sweep of funding increase rate (limited funding condition).
        """
        incr = 0.02                 # yearly increase start value
        incrby = 0.005              # increase increment for each loop

        aparms = self._load_data_file('alife_sweep3')
        lparms = {**self.parms, **aparms}

        # Parameter values for all simulation loops
        for simul in lparms['loops']:
            lparms['prefix'] = base_prefix + "inc0" + str(simul) + "_"
            lparms['yearly_increase'] = incr
            self._run_sim(dfile, lparms, lparms['prefix'])
            incr += incrby          # update yearly increase

        return lparms

    # ------------------------------------------------------------
    # Alife sweeps 4
    # ------------------------------------------------------------

    def run_alife_sweep4(self, dfile, base_prefix):
        """
        Run sweep of funding limits -- baseline scenario vs postdocs.
        """

        # run funding limited sim - baseline
        aparms = self._load_data_file('alife_sweep4')
        lparms1 = {**self.parms, **aparms['sim1']}
        print('Funding Limited -- Baseline')
        self._run_sim(dfile, lparms1, base_prefix)

        # run funding limited sim - postdocs
        lparms2 = {**lparms1, **aparms['sim2']}
        print('Funding Limited -- Postdocs')
        self._run_sim(dfile, lparms2, base_prefix)
        return lparms2

    # ------------------------------------------------------------
    # Alife sweeps 5
    # ------------------------------------------------------------

    def run_alife_sweep5(self, dfile, base_prefix):
        """
        Run sweep of funding limits -- baseline scenario vs postdocs.
        """

        indx = 0
        pdocs = 0
        aparms = self._load_data_file('alife_sweep5')
        lparms = {**self.parms, **aparms}

        for simul in [
                "limited_base", "limited_pdr", "unlimited_base",
                "unlimited_pdr"]:
            lparms['growing_pop'] = 1 - pdocs
            lparms['limited_funding'] = (indx < 2)
            lparms['prefix'] = base_prefix + "funding_" + simul
            lparms['use_postdocs'] = pdocs
            print(simul)
            self._run_sim(dfile, lparms, lparms['prefix'])
            indx += 1
            pdocs = 1 - pdocs

        return lparms

    # ------------------------------------------------------------
    # Alife thermostat runs
    # ------------------------------------------------------------

    def run_alife_thermo(self, dfile, base_prefix):
        """
        Run sims for Alife XV paper with thermostat.
        """
        indx = 0
        pdocs = 0
        aparms = self._load_data_file('alife_thermo')
        lparms = {**self.parms, **aparms}

        for simul in ["growingPop", "noRQnoM", "RQnoM", "RQM"]:
            lparms['growing_pop'] = 1 - pdocs
            lparms['mentored_pdrs'] = indx / 3
            lparms['pdr_rq_counts'] = (indx / 2) - (indx / 4)
            lparms['prefix'] = base_prefix + "thermo_" + simul
            lparms['use_postdocs'] = pdocs
            print(simul)
            self._run_sim(dfile, lparms, lparms['prefix'])
            indx += 1
            if indx == 1:
                pdocs = 1 - pdocs

        return lparms

    # ------------------------------------------------------------
    # Run gem sweep
    # ------------------------------------------------------------

    def run_gem_sweep1(self):
        """
        Runs for sensitivity analysis using GEM-SA
        """
        # Abbreviations
        # 'r_' prefix: research
        # '_f' suffix: file
        r_means = []
        r_ses = []

        gparms = self._load_data_file('gem_sweep1')
        lparms = {**self.parms, **gparms}

        dataf = open(self.prefs['roi_file'], 'w')
        gemf = open(self.prefs['gdata_file'], 'w')
        meansf = open(self.prefs['gmeans_file'], 'w')
        outf = open(self.prefs['gout_file'], 'w')
        prefix = lparms['prefix']
        sim_runs = lparms['runs']

        _headers(dataf, True, True, True, False)

        for achance in lparms['chance']:
            for amentor in lparms['bonus2']:
                for anewb in lparms['time1']:
                    for ahunt in lparms['time1']:
                        lparms['jobhunt_time'] = ahunt
                        lparms['mentoring_bonus'] = amentor
                        lparms['newb_time'] = anewb
                        lparms['postdoc_chance'] = achance
                        print("Trying postdoc chance: ", achance)
                        print("Trying mentoring bonus: ", amentor)
                        print("Trying newb time: ", anewb)
                        print("Trying jobhunting time: ", ahunt)
                        r_list = []
                        r_sum = 0.0
                        meansf.write(_tabs([achance, amentor, anewb, ahunt]))
                        for arun in range(0, sim_runs):
                            print(arun)
                            r_out = self._run_sim(
                                dataf, lparms, prefix, sim=True)[2]
                            r_list.append(r_out)
                            r_sum += r_out
                            print(r_out)
                            gemf.write(_tabs(
                                [achance, amentor, anewb, ahunt, r_out]))
                            r_means.append(pylab.mean(r_list))
                            outf.write(str(r_sum/sim_runs) + "\n")
                            r_ses.append(pylab.std(r_list) / \
                                         math.sqrt(sim_runs))

        dataf.close()
        gemf.close()
        meansf.close()
        outf.close()

    # ------------------------------------------------------------
    # Run gem sweep
    # ------------------------------------------------------------

    def run_gem_sweep2(self):
        """
        Runs for sensitivity analysis using GEM-SA -- NO MENTORING
        """
        # 'r_' prefix indicates 'research'
        r_means = []
        r_ses = []

        gparms = self._load_data_file('gem_sweep2')
        lparms = {**self.parms, **gparms}

        dataf = open(self.prefs['roi_file'], 'w')
        gemf = open(self.prefs['gdata_file'], 'w')
        meansf = open(self.prefs['gmeans_file'], 'w')
        outf = open(self.prefs['gout_file'], 'w')
        prefix = lparms['prefix']
        sim_runs = lparms['runs']

        _headers(dataf, True, True, True, False)

        for achance in lparms['chance']:
            for anewb in lparms['time2']:
                for ahunt in lparms['time2']:
                    lparms['jobhunt_time'] = ahunt
                    lparms['newb_time'] = anewb
                    lparms['postdoc_chance'] = achance
                    print("Trying postdoc chance: ", achance)
                    print("Trying newb time: ", anewb)
                    print("Trying jobhunting time: ", ahunt)
                    r_list = []
                    r_sum = 0.0
                    meansf.write(_tabs([achance, anewb, ahunt]))
                    for arun in range(0, sim_runs):
                        print(arun)
                        r_out = self._run_sim(
                            dataf, lparms, prefix, sim=True)[2]
                        r_list.append(r_out)
                        r_sum += r_out
                        print(r_out)
                        gemf.write(_tabs([achance, anewb, ahunt, r_out]))
                        r_means.append(pylab.mean(r_list))
                        outf.write(str(r_sum/sim_runs) + "\n")
                        r_ses.append(pylab.std(r_list) / \
                                     math.sqrt(sim_runs))

        dataf.close()
        gemf.close()
        meansf.close()
        outf.close()

    # ------------------------------------------------------------
    # Run gem sweep
    # ------------------------------------------------------------

    def run_gem_sweep3(self):
        """
        Runs for sensitivity analysis using GEM-SA --
        INCLUDING MENTORING SWITCH
        """
        r_means = []
        r_ses = []

        gparms = self._load_data_file('gem_sweep3')
        lparms = {**self.parms, **gparms}

        dataf = open(self.prefs['roi_file'], 'w')
        gemf = open(self.prefs['gem_file'], 'w')
        meansf = open(self.prefs['means_file'], 'w')
        outf = open(self.prefs['out_file'], 'w')
        prefix = lparms['prefix']
        sim_runs = lparms['runs']

        _headers(dataf, True, True, True, False)

        for achance in lparms['chance']:
            for amentor in lparms['bonus2']:
                for anewb in lparms['time1']:
                    for ahunt in lparms['time1']:
                        lparms['postdoc_chance'] = achance
                        lparms['mentoring_bonus'] = amentor
                        lparms['newb_time'] = anewb
                        lparms['jobhunt_time'] = ahunt
                        print("Trying postdoc chance: ", achance)
                        print("Trying mentoring bonus: ", amentor)
                        print("Trying newb time: ", anewb)
                        print("Trying jobhunting time: ", ahunt)
                        r_list = []
                        r_sum = 0.0
                        meansf.write(_tabs(
                            [achance, amentor, anewb, ahunt]))
                        for arun in range(0, sim_runs):
                            print(arun,)
                            r_out = self._run_sim(
                                dataf, lparms, prefix, sim=True)[2]
                            r_list.append(r_out)
                            r_sum += r_out
                            print(r_out)
                            gemf.write(_tabs(
                                [achance, amentor, anewb, ahunt, r_out]))
                            r_means.append(pylab.mean(r_list))
                            outf.write(str(r_sum/sim_runs) + "\n")
                            r_ses.append(
                                pylab.std(r_list) / math.sqrt(sim_runs))

        dataf.close()
        gemf.close()
        meansf.close()
        outf.close()

    # ------------------------------------------------------------
    # Run gem sweep
    # ------------------------------------------------------------

    def run_gem_sweep4(self):
        """
        Runs for sensitivity analysis using GEM-SA --
        INCLUDING MENTORING SWITCH
        """
        r_means = []
        r_ses = []

        gparms = self._load_data_file('gem_sweep4')
        lparms = {**self.parms, **gparms}

        dataf = open(self.prefs['roi_file'], 'w')
        gemf = open(self.prefs['gem_file'], 'w')
        meansf = open(self.prefs['means_file'], 'w')
        outf = open(self.prefs['out_file'], 'w')
        prefix = lparms['prefix']
        sim_runs = lparms['runs']

        _headers(dataf, True, True, True, False)

        for achance in lparms['chance']:
            for amentor in lparms['bonus2']:
                for anewb in lparms['time1']:
                    for ahunt in lparms['time1']:
                        lparms['postdoc_chance'] = achance
                        lparms['mentoring_bonus'] = amentor
                        lparms['newb_time'] = anewb
                        lparms['jobhunt_time'] = ahunt
                        print("Trying postdoc chance: ", achance)
                        print("Trying mentoring bonus: ", amentor)
                        print("Trying newb time: ", anewb)
                        print("Trying jobhunting time: ", ahunt)
                        redundl = []
                        redunds = 0.0
                        meansf.write(_tabs(
                            [achance, amentor, anewb, ahunt]))
                        for arun in range(0, sim_runs):
                            print(arun,)
                            #sim = self._run_sim(
                            #dataf, lparms, prefix, sim=True)
                            redundo = self._run_sim(
                                dataf, lparms, prefix, sim=True)[3]
                            redundl.append(redundo)
                            redunds += redundo
                            print(redundo)
                            gemf.write(_tabs(
                                [achance, amentor, anewb, ahunt, redundo]))
                            r_means.append(pylab.mean(redundl))
                            outf.write(str(redunds/sim_runs) + "\n")
                            r_ses.append(
                                pylab.std(redundl) / math.sqrt(sim_runs))

        dataf.close()
        gemf.close()
        meansf.close()
        outf.close()

    # ------------------------------------------------------------
    # Run gem sweep
    # ------------------------------------------------------------

    def run_gem_sweep5(self):
        """
        Runs for sensitivity analysis using GEM-SA --
        INCLUDING MENTORING SWITCH
        """
        r_means = []
        r_ses = []

        gparms = self._load_data_file('gem_sweep5')
        lparms = {**self.parms, **gparms}

        dataf = open(self.prefs['roi_file'], 'w')
        gemf = open(self.prefs['gem_file'], 'w')
        meansf = open(self.prefs['means_file'], 'w')
        outf = open(self.prefs['out_file'], 'w')
        prefix = lparms['prefix']
        sim_runs = lparms['runs']

        _headers(dataf, True, True, True, False)

        for achance in lparms['chance']:
            for amentor in lparms['bonus1']:
                for astartgrant in lparms['grant']:
                    for aninc in lparms['increase']:
                        lparms['postdoc_chance'] = achance
                        lparms['mentoring_bonus'] = amentor
                        lparms['starting_grant_fund'] = astartgrant
                        lparms['yearly_increase'] = aninc
                        redundl = []
                        redunds = 0.0
                        meansf.write(
                            _tabs([achance, amentor, astartgrant, aninc]))
                        for arun in range(0, sim_runs):
                            print(arun,)
                            #sim = self._run_sim(prefix, sim=True)
                            redundo = self._run_sim(
                                dataf, lparms, prefix, sim=True)[3]
                            redundl.append(redundo)
                            redunds += redundo
                            print(redundo)
                            gemf.write(_tabs(
                                [achance, amentor, astartgrant, aninc,
                                 redundo]))
                            r_means.append(pylab.mean(redundl))
                            outf.write(str(redunds/sim_runs) + "\n")
                            r_ses.append(
                                pylab.std(redundl) / math.sqrt(sim_runs))

        dataf.close()
        gemf.close()
        meansf.close()
        outf.close()

    # ------------------------------------------------------------
    # Run gem sweep
    # ------------------------------------------------------------

    def run_gem_sweep6(self):
        """
        Runs for sensitivity analysis using GEM-SA --
        INCLUDING limited funding switch
        """
        r_means = []
        r_ses = []

        gparms = self._load_data_file('gem_sweep5')
        lparms = {**self.parms, **gparms}

        dataf = open(self.prefs['roi_file'], 'w')
        gemf = open(self.prefs['gem_file'], 'w')
        meansf = open(self.prefs['means_file'], 'w')
        outf = open(self.prefs['out_file'], 'w')
        prefix = lparms['prefix']
        sim_runs = lparms['runs']
        limitfunding = [True, False]

        _headers(dataf, True, True, True, False)

        for achance in lparms['chance']:
            for amentor in lparms['bonus1']:
                for ahunt in lparms['time2']:
                    for afund in limitfunding:
                        lparms['postdoc_chance'] = achance
                        lparms['mentoring_bonus'] = amentor
                        lparms['jobhunt_time'] = ahunt
                        lparms['limited_funding'] = afund
                        redundl = []
                        redunds = 0.0
                        meansf.write(_tabs(
                            [achance, amentor, ahunt, afund]))
                        for arun in range(0, sim_runs):
                            print(arun,)
                            # sim = self._run_sim(
                            # dataf, lparms, prefix, sim=True)
                            redundo = self._run_sim(
                                dataf, lparms, prefix, sim=True)[3]
                            redundl.append(redundo)
                            redunds += redundo
                            print(redundo)
                            gemf.write(_tabs(
                                [achance, amentor, ahunt, afund, redundo]))
                            r_means.append(pylab.mean(redundl))
                            outf.write(str(redunds/sim_runs) + "\n")
                            r_ses.append(
                                pylab.std(redundl) / math.sqrt(sim_runs))

        dataf.close()
        gemf.close()
        meansf.close()
        outf.close()

    # ------------------------------------------------------------
    # Run a single run
    # ------------------------------------------------------------

    def run_once(self):
        """
        Run a single simulation
        Note no use of lparms as used as read-only
        """
        dfile = open(self.prefs['roi_file'], 'w')
        _headers(dfile, False, False, False, False)
        self._run_sim(
            dfile, self.parms, self.parms['prefix'], sim=True)
        dfile.close()

    # ------------------------------------------------------------
    # Run partial *** needs to be looked at ***
    # ------------------------------------------------------------

    # Disable groupby checking in Pylint
    # pylint: disable=E1101

    def run_partial(self):
        """
        Run partial -- needs to be looked at for functionality
        """
        plt.clf()
        fields = self.parms['pdrgrow'] + self.parms['prqpdrs'] + \
                 self.parms['totrosack']
        roif = pd.read_csv(self.prefs['roi_file'], usecols=fields)
        deviation = roif.groupby(
            self.parms['pdrgrow'] + self.parms['prqpdrs']).std()
        means = roif.groupby(
            self.parms['pdrgrow'] + self.parms['prqpdrs']).mean()
        print(means)
        print(deviation)
        scenariol = ['Growing Pop', 'NoRQnoM', 'noRQM', 'RQnoM', 'RQM']

        # Total RO plot
        fig = _plot_fig(
            False, means[TOTRO], 'Output in Postdoc Scenarios',
            deviation.Total_RO, self.prefs['palettes'][0], "Scenarios",
            "Research Output", scenariol, self.parms['prefix'],
            MEMRUN_TOT, False)

        # Mean Output plot
        _plot_fig(
            fig, means[MEANROALL], 'Mean Output in Postdoc Scenarios',
            deviation.Mean_RO_all, self.prefs['palettes'][0],
            "Scenarios", "Mean Research Output",
            scenariol, self.parms['prefix'], MEMRUN_MEAN)

        # ROI in Postdoc Scenarios plot
        _plot_fig(
            fig, means['ROI_noPDRs'], 'ROI in Postdoc Scenarios',
            deviation.ROI_no_PDRs, 'g', "Scenarios",
            "Return on Investment",
            scenariol, self.parms['prefix'], MEMRUN_ROI)

        datafr = roif[roif.growing_pop != 1]
        deviationr = datafr.groupby(
            [self.parms['usepdrs'], self.parms['prqpdrs']]).std()
        meansr = datafr.groupby(
            [self.parms['usepdrs'], self.parms['prqpdrs']]).mean()
        scenariolr = ['NoRQnoM', 'noRQM', 'RQnoM', 'RQM']

        # Redundancies in Postdoc Scenarios plot
        _plot_fig(
            fig, meansr[TOTSACK], 'Redundancies in Postdoc Scenarios',
            deviationr.Total_Sacked, 'r', "Scenarios", "Redundancies",
            scenariolr, self.parms['prefix'], MEMRUN_SACK)


    # ------------------------------------------------------------
    # Run sweep
    # ------------------------------------------------------------

    def run_sweep1(self):
        """
        Description of sweep
        """
        dataf = open(self.prefs['roi_file'], 'w')
        _headers(dataf, False, True, True, True)
        prefix = self.parms['prefix']
        self.run_alife_sweep1(dataf, self.parms['prefix'])
        dataf.close()
        plt.clf()
        fields = self.parms['fundlim'] + self.parms['totrosack']
        fundinglist = [True, False]

        # Fill data frame from CSV file
        roif = pd.read_csv(self.prefs['roi_file'], usecols=fields)
        deviation = roif.groupby(FUNDLIM).std()
        means = roif.groupby(FUNDLIM).mean()
        print(means)
        print(deviation)

        # Limited Funding vs Output plot
        fig = _plot_fig(
            False, means[TOTRO], 'Limited Funding vs Output',
            deviation.Total_RO, 'c', "Funding Limit", "Research Output",
            fundinglist, prefix, PRORUN_TOT, False)

        # Limited Funding vs Mean Output plot
        _plot_fig(
            fig, means['Mean_R0_all'], 'Limited Funding vs Mean Output',
            deviation.Mean_RO_all, 'm', "Funding Limit",
            "Mean Research Output",
            fundinglist, prefix, PRORUN_MEAN)

        # Limited Funding vs ROI plot
        _plot_fig(
            fig, means[ROINOPDRS], 'Limited Funding vs ROI',
            deviation.ROI_no_PDRs, 'g', "Funding Limit",
            "Return on Investment",
            fundinglist, prefix, PRORUN_ROI)

        plt.clf()
        #datafr = dataf[roif.promo_chance != 1]
        #meansr = datafr.groupby(self.parms['prochan']).mean()
        #deviationr = datafr.groupby(self.parms['prochan']).std()

        # Limited Funding vs Redundancies plot
        _plot_fig(
            fig, means[TOTSACK], 'Limited Funding vs Redundancies',
            deviation.Total_Sacked, 'r', "Funding Limit",
            "Redundancies", fundinglist, prefix, PRORUN_SACK, False)


    # ------------------------------------------------------------
    # Run sweep
    # ------------------------------------------------------------

    def run_sweep2(self):
        """
        Description of sweep
        """
        # Alife XV promo chance sweep
        dataf = open(self.prefs['roi_file'], 'w')
        _headers(dataf, False, True, True, False)
        prefix = self.parms['prefix']
        self.run_alife_sweep2(dataf, prefix)
        dataf.close()
        plt.clf()

        fields = self.parms['prochan'] + self.parms['totrosack']
        roif = pd.read_csv(self.prefs['roi_file'], usecols=fields)
        deviation = roif.groupby(PROCHAN).std()
        means = roif.groupby(PROCHAN).mean()
        promochancelist = [15, 25, 50, 75, 100]
        print(means)
        print(deviation)

        # Promotion chance vs Output plot
        fig = _plot_fig(
            False, means[TOTRO], 'Promotion Chance vs Output',
            deviation.Total_RO, 'c', "Promotion Chance",
            "Research Output",
            promochancelist, prefix, PRORUN_TOT)

        # Promotion chance vs Mean Output plot
        _plot_fig(
            fig, means[MEANROALL], 'Promotion Chance vs Mean Output',
            deviation.Mean_RO_all, 'm', "Promotion Chance",
            "Mean Research Output",
            promochancelist, prefix, PRORUN_MEAN)

        # Promotion chance vs Mean Output plot
        _plot_fig(
            fig, means[ROINOPDRS], 'Promotion Chance vs ROI',
            deviation.ROI_no_PDRs, 'g', "Promotion Chance",
            "Return on Investment",
            promochancelist, prefix, PRORUN_ROI)

        datafr = dataf[roif.promo_chance != 1]
        deviationr = datafr.groupby(self.parms['prochan']).std()
        meansr = datafr.groupby(self.parms['prochan']).mean()

        # Redundancies vs Promotion chance plot
        _plot_fig(
            fig, meansr[TOTSACK], 'Redundancies vs Promotion Chance',
            deviationr.Total_Sacked, 'r', "Promotion Chance",
            "Redundancies",
            promochancelist, prefix, PRORUN_SACK)


    # ------------------------------------------------------------
    # Run sweep
    # ------------------------------------------------------------

    def run_sweep3(self):
        """
        Description of sweep
        """
        dataf = open(self.prefs['roi_file'], 'w')
        _headers(dataf, False, True, True, True)
        prefix = self.parms['prefix']
        self.run_alife_sweep3(dataf, prefix)
        dataf.close()
        plt.clf()
        fields = self.parms['fundinc'] + self.parms['totrosack']
        roif = pd.read_csv(self.prefs['roi_file'], usecols=fields)
        deviation = roif.groupby(FUNDINC).std()
        fundinglist = [2, 2.5, 3, 3.5, 4, 4.5, 5]
        means = roif.groupby(FUNDINC).mean()
        print(means)
        print(deviation)

        # Funding Increases vs Output plot
        fig = _plot_fig(
            False, means[TOTRO], 'Funding Increases vs Output',
            deviation.Total_RO, 'c', "Funding Increase",
            "Research Output",
            fundinglist, prefix, INCRUN_TOT, False)

        # Funding Increases vs Mean Output
        _plot_fig(
            fig, means[MEANROALL], 'Funding Increases vs Mean Output',
            deviation.Mean_RO_all, 'm', "Funding Increase",
            "Mean Research Output",
            fundinglist, prefix, INCRUN_MEAN)

        # Funding Increases vs ROI
        _plot_fig(
            fig, means[ROINOPDRS], 'Funding Increases vs ROI',
            deviation.ROI_no_PDRs, 'g', "Funding Increase",
            "Return on Investment",
            fundinglist, prefix, INCRUN_ROI)

        #datafr = dataf[roif.promo_chance != 1]
        #meansr = datafr.groupby(self.parms['prochan']).mean()
        #deviationr = datafr.groupby(self.parms['prochan']).std()

        # Funding Increases vs Redundancies
        _plot_fig(
            fig, means[TOTSACK], 'Funding Increases vs Redundancies',
            deviation.Total_Sacked, 'r', "Funding Increase",
            "Redundancies",
            fundinglist, prefix, INCRUN_SACK)

    # ------------------------------------------------------------
    # Run sweep
    # ------------------------------------------------------------

    def run_sweep4(self):
        """
        Description of sweep
        """
        # Parameter sweep, baseline vs limited vs unlimited funding
        dataf = open(self.prefs['roi_file'], 'w')
        _headers(dataf, False, True, True, True)
        prefix = self.parms['prefix']
        self.run_alife_sweep4(dataf, prefix)
        dataf.close()
        plt.clf()
        fields = self.parms['pdrgrow'] + self.parms['totrosack']
        roif = pd.read_csv(self.prefs['roi_file'], usecols=fields)
        deviation = roif.groupby(self.parms['pdrgrow']).std()
        means = roif.groupby(self.parms['pdrgrow']).mean()
        print(means)
        print(deviation)
        scenariol = ['Growing Pop', 'Limited funding']

        # Output Comparison: Baseline vs Limited Funding
        fig = _plot_fig(
            False, means[TOTRO],
            'Output Comparison: Baseline vs Limited Funding',
            deviation.Total_RO, self.prefs['palettes'][1], "Scenarios",
            "Research Output", scenariol, prefix, MEMRUN_TOT, False)

        # Mean Output Comparison: Baseline vs Limited Funding
        _plot_fig(
            fig, means[TOTSACK],
            'Mean Output Comparison: Baseline vs Limited Funding',
            deviation.Mean_RO_all, self.prefs['palettes'][2],
            "Scenarios",
            "Main Research Output", scenariol, prefix, MEMRUN_MEAN)

        # ROI in Limited-Funding Postdoc Scenarios
        _plot_fig(
            fig, means[ROINOPDRS],
            'ROI in Limited-Funding Postdoc Scenarios',
            deviation.ROI_no_PDRs, 'g', "Scenarios",
            "Return on Investment",
            scenariol, prefix, MEMRUN_ROI)

        # Redundancies: Baseline vs Limited Funding
        _plot_fig(
            fig, means[TOTSACK],
            'Redundancies: Baseline vs Limited Funding',
            deviation.ROI_no_PDRs, 'r', "Scenarios", "Redundancies",
            scenariol, prefix, INCRUN_SACK)

    # ------------------------------------------------------------
    # Run sweep
    # ------------------------------------------------------------

    def run_sweep5(self):
        """
        Description of sweep
        """
        # Parameter sweep, baseline vs limited vs unlimited funding
        dataf = open(self.prefs['roi_file'], 'w')
        _headers(dataf, False, True, True, True)
        prefix = self.parms['prefix']
        self.run_alife_sweep5(dataf, prefix)
        dataf.close()
        plt.clf()
        fields = self.parms['pdrgrow'] + self.parms[
            'totrosack'] + self.parms['fundlim']
        roif = pd.read_csv(self.prefs['roi_file'], usecols=fields)
        deviation = roif.groupby(
            self.parms['pdrgrow'] + self.parms['fundlim']).std()
        means = roif.groupby(
            self.parms['pdrgrow'] + self.parms['fundlim']).mean()
        print(means)
        print(deviation)
        scenariol = ['GP-LF', 'PDR-LF', 'GP-UF', 'PDR-UF']

        # Output Comparison: Baseline vs Limited Funding
        fig = _plot_fig(
            False, means[TOTRO],
            'Output Comparison: Baseline vs Limited Funding',
            deviation.Total_RO, self.prefs['palettes'][3], "Scenarios",
            "Research Output", scenariol, prefix, MEMRUN_TOT, False)

        # Mean Output Comparison: Baseline and Funding Scenarios
        _plot_fig(
            fig, means[TOTSACK],
            'Mean Output Comparison: Baseline and Funding Scenarios',
            deviation.Mean_RO_all, self.prefs['palettes'][3],
            "Scenarios",
            "Mean Researhc Output", scenariol, prefix, MEMRUN_MEAN)

        # Mean Output Comparison: Baseline and Funding Scenarios
        _plot_fig(
            fig, means[ROINOPDRS], 'ROI by Scenario',
            deviation.ROI_nno_PDRs, 'g', "Scenarios",
            "Mean Research Output",
            scenariol, prefix, MEMRUN_MEAN)

        # ROI by Scenario
        _plot_fig(
            fig, means[ROINOPDRS], 'ROI by Scenario',
            deviation.ROI_nno_PDRs, 'g', "Scenarios",
            "Return on Investment",
            scenariol, prefix, MEMRUN_ROI)

        # Redundancies: Baseline and Funding Scenarios
        _plot_fig(
            fig, means[TOTSACK],
            'Redundancies: Baseline and Funding Scenarios',
            deviation.Total_Sacked, 'r', "Scenarios", "Redundancies",
            scenariol, prefix, MEMRUN_SACK)

        # pylint: enable=E1101

    # ------------------------------------------------------------
    # WCSS runs
    # ------------------------------------------------------------

    def run_wcss_sims(self, dfile):
        """
        Run basic simulations associated with the WCSS paper.
        """
        wparms = self._load_data_file('wcss_sims')
        lparms1 = {**self.parms, **wparms['sim1']}
        print('thermostat')
        self._run_sim(dfile, lparms1, lparms1['prefix'])

        # run MEMORY A model ("bad" parameters)
        lparms2 = {**lparms1, **wparms['sim2']}
        print('memory A')
        self._run_sim(dfile, lparms2, lparms2['prefix'])

        # run MEMORY B model ("good" parameters)
        lparms3 = {**lparms2, **wparms['sim3']}
        print('memory B')
        self._run_sim(dfile, lparms3, lparms3['prefix'])

        # run FIXED model
        lparms4 = {**lparms3, **wparms['sim4']}
        print('fixed')
        self._run_sim(
            dfile, lparms4, lparms4['prefix'], sim=False)


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
