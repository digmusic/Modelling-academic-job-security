"""
--------------------------------------------------------------------------------
FILE: Dispatch.py
DESC: Dispatch simulations
AUTH: thorsilver
VERS: 1.0 March 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""

# Unused imports *** check that this is correct ***
# from collections import defaultdict
# from random import Random
# import csv
# import os
# import numpy

from __future__ import print_function # ensure Python 3.x print form known
import math
import Utils
import matplotlib.pyplot as plt
import pandas as pd
import pylab
from Simulation import Simulation

DIR_PREFIX = 'run-files/'      # subdir of current dir *** testing only ***
# Original value was as follows:
# '/Users/u0030612/Documents/results/test/Alife/limited funding/baseComp8/'

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

# Filename constants with omitted extensions (usually .csv or .txt files)
ROI_TEST = "roi_test"
SUMM_DATA = "summ_data"
SUMM_STATS = "summ_stats"
GEMSA_DATA = "GEMSA data new"
GEMSA_MEANS = "GEMSA means new"
GEMSA_OUTS = "GEMSA outputs new"

# Field names: level 0
FUNDINC = 'Funding_Increase'
FUNDLIM = 'Funding_Limit'
GROWPOP = 'Growing_Pop'
MEANROALL = 'Mean_RO_all'
PROCHAN = 'Promo_Chance'
ROINOPDRS = 'ROI_no_PDRs'
TOTRO = 'Total_RO'
TOTSACK = 'Total_Sacked'
USEPDRS = 'Use_PDRs'

# Field names: level 1
FUNDINC_L = [FUNDINC]
FUNDLIM_L = [FUNDLIM]
GROWPOP_L = [GROWPOP]
PROCHAN_L = [PROCHAN]
RQPDRS_L = ['RQ_counts', 'Mentored_PDRs']
TOTROSACK_L = [TOTRO, MEANROALL, ROINOPDRS, TOTSACK]
USEPDRS_L = [USEPDRS]

# Field names: level 2
PDRGROW_L = USEPDRS_L + GROWPOP_L

# Colour palettes
PALETTE1_L = ['Orange', 'Salmon', 'Violet', 'Chocolate', 'LightSteelBlue']
PALETTE2_L = ['Orange', 'Salmon']
PALETTE3_L = ['Orange', 'Salmon', 'Violet']
PALETTE4_L = ['Orange', 'Salmon', 'Violet', 'Green']

# Loop values for runs
BONUS1_L = [0.30, 0.40, 0.50, 0.60]
BONUS2_L = [0.30, 0.40, 0.50, 0.60, 0.70]
CHANCE_L = [0.15, 0.25, 0.50, 0.75, 1.0]
GRANT_L = [10, 20, 30, 40, 50]
INCREASE_L = [0.01, 0.02, 0.03, 0.04]
TIME1_L = [0.1, 0.3, 0.5, 0.7]
TIME2_L = [0.1, 0.3, 0.5, 0.7, 0.9]


# ============================================================
# Private helper functions to remove (largely) duplicate code
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


def _plot_fig(depfig, item, ptitle, errval, fcolor, xlabel, ylabel,
              ticklbl, rootdir, fname, clflag=True):

    """
    Plot a figure
    """
    if clflag:
        plt.clf()
    fig1 = item.plot(
        kind='bar', title=ptitle, yerr=errval, facecolor=fcolor, alpha=0.5, rot=0)
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


def _run_sim(dfile, parms, exp_prefix, rng, sim=True):

    """
    FUNC: Run multiple initialisations of single set of experimental conditions.
    ARGS: dfile      - data file
          parms      - parameters to control the runs
          exp_prefix - prefix for experiments
          rng        - range
          sim        - simulation flag. If sim=False, output is estimated,
                       rather than simulated.  Currently, this is just used for
                       evaluating the "fixed time investment" heuristic.
    RETS:
    """
    sum_roi_values = []
    roi_values = []
    total_r_values = []
    mean_tg_values = []
    corr_rq_held_values = []
    redundancies_final = 0

    # for each initialisation
    for run in range(parms['runs']):
        parms['seed'] = rng.randint(0, 99999999)
        parms['prefix'] = exp_prefix + "%d/" % parms['seed']
        print("   run #%d" % run)

        # run experiment; either simulation or estimation
        exp = Simulation(parms)
        if sim:
            exp.run()
        else:
            exp.test_flat(parms['fixed_time'])

        roi_values.append(exp.calc_roi(1)[2])
        total_r_values.append(exp.calc_mean_total_output())
        mean_tg_values.append(exp.calc_mean_time_grant())
        corr_rq_held_values.append(exp.calc_mean_corr_rq_held())
        sum_roi_values.append(exp.calc_roi_sum)
        #roi_SEs.append(numpy.std(sum_roi_values))
        exp.calc_redundancies()
        redundancies_final = exp.redundancies_total

        if parms['write_output']:
            Utils.create_dir(parms['prefix'])
            exp.write_output()
            #redundancies_final = exp.redundancies_total
            #redundancies_final.append(exp.redundancies_total)
            dfile.write(
                str(run + 1) + "," + str(parms['use_postdocs']) + "," + \
                str(parms['growing_pop']) + "," + \
                str(parms['pdr_rq_counts']) + "," + \
                str(parms['mentored_pdrs']) + "," + \
                str(exp.roi_sum[parms['iterations'] - 1]) + "," + \
                str(exp.roi_sum_pdr[parms['iterations'] - 1]) + "," + \
                str(total_r_values[run]) + "," + \
                str(exp.mean_r[parms['iterations'] - 1]) + "," + \
                str(exp.mean_r_old_academic[parms['iterations'] - 1]) + "," + \
                str(exp.mean_r_postdoc[parms['iterations'] - 1]) + "," + \
                str(exp.mean_r_former_pdr[parms['iterations'] - 1]) + "," + \
                str(mean_tg_values[run]) + "," + \
                str(parms['learning_type']) + "," + \
                str(parms['postdoc_chance']) + "," + \
                str(exp.redundancies_total) + "," + \
                str(parms['limited_funding']) + "," + \
                str(parms['yearly_increase']) + "\n")

    # calculate and write summary statistics
    summ_data = [roi_values, total_r_values, mean_tg_values,
                 corr_rq_held_values]

    summ_stats = [pylab.mean(roi_values), pylab.std(roi_values),
                  pylab.mean(total_r_values), pylab.std(total_r_values),
                  pylab.mean(mean_tg_values), pylab.std(mean_tg_values),
                  pylab.mean(corr_rq_held_values),
                  pylab.std(corr_rq_held_values)]

    if parms['write_output']:
        Utils.write_data_2d(summ_data, exp_prefix + SUMM_DATA + '.csv')
        Utils.write_data(summ_stats, exp_prefix + SUMM_STATS + '.csv')

    return pylab.mean(roi_values), pylab.std(roi_values), \
        total_r_values[-1], redundancies_final


def _tab_str3n(str1, str2, str3):

    """
    Construct a string of the arguments, separated by tabs and
    finished with a newline
    """
    return str(str1) + "\t" + str(str2) + "\t" + str(str3) + "\n"

def _tab_str4(str1, str2, str3, str4):

    """
    Construct a string of the arguments, separated by tabs and
    finished with a newline
    """
    return str(str1) + "\t" + str(str2) + "\t" + str(str3) + "\t" + str(str4)

def _tab_str4n(str1, str2, str3, str4):

    """
    Construct a string of the arguments, separated by tabs and
    finished with a newline
    """
    return _tab_str4(str1, str2, str3, str4) + "\n"


def _tab_str5n(str1, str2, str3, str4, str5):

    """
    Construct a string of the arguments, separated by tabs and
    finished with a newline
    """
    return _tab_str4(str1, str2, str3, str4) + str(str5) + "\n"


# ============================================================
# Public functions: level 0
# ============================================================

def alife_run(func, parms, seedrng, fund):

    """
    Run an Alife XV-based simulation.
    Expects the function to be called using parms, parms-prefix, and seedrng
    """
    dfile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(dfile, False, True, True, fund)
    rdir = parms['prefix']           # root directory
    func(parms, parms['prefix'], seedrng)
    dfile.close()
    return rdir

# ============================================================
# Public functions: level 1
# ============================================================

# ------------------------------------------------------------
# Alife sim runs
# ------------------------------------------------------------

# *** The code extracts parms that do not appear to change after initial value is
#     set and places them outside of the loop. If this assumption is incorrect
#     then the sims functions need to be looked at
# ***

def run_alife_sims(dfile, parms, base_prefix, seedrng):
    """
    Run sims for Alife XV paper.
    """
    parms['learning_type'] = 'memory'
    parms['runs'] = 50
    pdocs = 0
    indx = 0
    for simul in ["growingPop", "noRQnoM", "RQnoM", "RQM", "noRQM"]:
        parms['growing_pop'] = 1 - pdocs
        parms['mentored_pdrs'] = indx / 3
        parms['pdr_rq_counts'] = (indx / 2) - (indx / 4)
        parms['prefix'] = base_prefix + "memory_" + simul
        parms['use_postdocs'] = pdocs
        print(simul)
        _run_sim(dfile, parms, parms['prefix'], seedrng)
        indx += 1
        if indx == 1:
            pdocs = 1 - pdocs
    return parms

# ------------------------------------------------------------
# Alife sweeps 1
# ------------------------------------------------------------

def run_alife_sweep1(dfile, parms, base_prefix, seedrng):
    """
    Run sweep of promotion chance values.
    """
    parms['growing_pop'] = 0
    parms['learning_type'] = 'memory'
    parms['mentored_pdrs'] = 0
    parms['pdr_rq_counts'] = 1
    parms['runs'] = 50
    parms['use_postdocs'] = 1
    for simul in [0.15, 0.25, 0.5, 0.75, 1]:
        pcent = int(100 * simul)
        parms['postdoc_chance'] = simul
        parms['prefix'] = base_prefix + "memory_chance" + pcent + "_"
        print('Memory, postdocs promo ' + pcent + 'percent')
        _run_sim(dfile, parms, parms['prefix'], seedrng)
    return parms

# ------------------------------------------------------------
# Alife sweeps 2
# ------------------------------------------------------------

def run_alife_sweep2(dfile, parms, base_prefix, seedrng):
    """
    Run sweep of funding limits.
    """
    parms['growing_pop'] = 0
    parms['learning_type'] = 'memory'
    parms['limited_funding'] = True
    parms['mentored_pdrs'] = 1
    parms['pdr_rq_counts'] = 1
    parms['postdoc_chance'] = 0.15
    parms['prefix'] = base_prefix + "funding_limited_"
    parms['runs'] = 100
    parms['use_postdocs'] = 1
    print('Funding Limited')
    _run_sim(dfile, parms, parms['prefix'], seedrng)

    # run funding unlimited sim
    parms['limited_funding'] = False
    parms['prefix'] = base_prefix + "funding_unlimited_"
    print('Funding Unlimited')
    _run_sim(dfile, parms, parms['prefix'], seedrng)
    return parms

# ------------------------------------------------------------
# Alife sweeps 3
# ------------------------------------------------------------

def run_alife_sweep3(dfile, parms, base_prefix, seedrng):
    """
    Run sweep of funding increase rate (limited funding condition).
    """
    incr = 0.02                 # yearly increase start value
    incrby = 0.005              # increase increment for each loop

    # Parameter values for all simulation loops
    parms['growing_pop'] = 0
    parms['limited_funding'] = True
    parms['mentored_pdrs'] = 0
    parms['pdr_rq_counts'] = 1
    parms['runs'] = 20
    parms['starting_grant_fund'] = 30
    parms['use_postdocs'] = 1

    for simul in [2, 25, 3, 35, 4, 45, 5]:
        parms['prefix'] = base_prefix + "inc0" + simul + "_"
        parms['yearly_increase'] = incr
        _run_sim(dfile, parms, parms['prefix'], seedrng)
        incr += incrby          # update yearly increase

    return parms

# ------------------------------------------------------------
# Alife sweeps 4
# ------------------------------------------------------------

def run_alife_sweep4(dfile, parms, base_prefix, seedrng):
    """
    Run sweep of funding limits -- baseline scenario vs postdocs.
    """

    # run funding limited sim - baseline
    parms['growing_pop'] = 1
    parms['learning_type'] = 'memory'
    parms['limited_funding'] = True
    parms['mentored_pdrs'] = 1
    parms['pdr_rq_counts'] = 1
    parms['postdoc_chance'] = 0.15
    parms['prefix'] = base_prefix + "funding_limited_base"
    parms['runs'] = 10
    parms['use_postdocs'] = 0
    print('Funding Limited -- Baseline')
    _run_sim(dfile, parms, parms['prefix'], seedrng)

    # run funding limited sim - postdocs
    parms['growing_pop'] = 0
    parms['prefix'] = base_prefix + "funding_limited_pdr"
    parms['use_postdocs'] = 1
    print('Funding Limited -- Postdocs')
    _run_sim(dfile, parms, parms['prefix'], seedrng)
    return parms

# ------------------------------------------------------------
# Alife sweeps 5
# ------------------------------------------------------------

def run_alife_sweep5(dfile, parms, base_prefix, seedrng):
    """
    Run sweep of funding limits -- baseline scenario vs postdocs.
    """

    parms['learning_type'] = 'memory'
    parms['mentored_pdrs'] = 1
    parms['pdr_rq_counts'] = 1
    parms['postdoc_chance'] = 0.15
    parms['runs'] = 50
    indx = 0
    pdocs = 0
    for simul in [
            "limited_base", "limited_pdr", "unlimited_base", "unlimited_pdr"]:
        parms['growing_pop'] = 1 - pdocs
        parms['limited_funding'] = (indx < 2)
        parms['prefix'] = base_prefix + "funding_" + simul
        parms['use_postdocs'] = pdocs
        print(simul)
        _run_sim(dfile, parms, parms['prefix'], seedrng)
        indx += 1
        pdocs = 1 - pdocs
    return parms

# ------------------------------------------------------------
# Alife thermostat runs
# ------------------------------------------------------------

def run_alife_thermo(dfile, parms, base_prefix, seedrng):
    """
    Run sims for Alife XV paper with thermostat.
    """
    parms['learning_type'] = 'thermostat'
    parms['runs'] = 25
    indx = 0
    pdocs = 0
    for simul in ["growingPop", "noRQnoM", "RQnoM", "RQM"]:
        parms['growing_pop'] = 1 - pdocs
        parms['mentored_pdrs'] = indx / 3
        parms['pdr_rq_counts'] = (indx / 2) - (indx / 4)
        parms['prefix'] = base_prefix + "thermo_" + simul
        parms['use_postdocs'] = pdocs
        print(simul)
        _run_sim(dfile, parms, parms['prefix'], seedrng)
        indx += 1
        if indx == 1:
            pdocs = 1 - pdocs

    return parms

# ------------------------------------------------------------
# Run gem sweep
# ------------------------------------------------------------

def run_gem_sweep1(parms, seedrng):

    """
    Runs for sensitivity analysis using GEM-SA
    """
    researchmeans = []
    researchses = []

    parms['write_output'] = False

    sim_runs = 20
    gemfile = open(parms['prefix'] + GEMSA_DATA + ".txt", 'w')
    meansfile = open(parms['prefix'] + GEMSA_MEANS + ".txt", 'w')
    outfile = open(parms['prefix'] + GEMSA_OUTS + ".txt", 'w')

    postdocchancelist = CHANCE_L
    mentoringbonuslist = BONUS2_L
    newbtimelist = TIME1_L
    jobhunttimelist = TIME1_L

    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, True, True, True, False)
    root_dir = parms['prefix']

    for achance in postdocchancelist:
        for amentor in mentoringbonuslist:
            for anewb in newbtimelist:
                for ajobhunt in jobhunttimelist:
                    parms['jobhunt_time'] = ajobhunt
                    parms['mentoring_bonus'] = amentor
                    parms['newb_time'] = anewb
                    parms['postdoc_chance'] = achance
                    print("Trying postdoc chance: ", achance)
                    print("Trying mentoring bonus: ", amentor)
                    print("Trying newb time: ", anewb)
                    print("Trying jobhunting time: ", ajobhunt)
                    researchlist = []
                    researchsum = 0.0
                    meansfile.write(_tab_str4n(achance, amentor, anewb, ajobhunt))
                    for arun in range(0, sim_runs):
                        print(arun)
                        researchout = _run_sim(datafile, parms, root_dir,
                                               seedrng, sim=True)[2]
                        researchlist.append(researchout)
                        researchsum += researchout
                        print(researchout)
                        gemfile.write(
                            _tab_str5n(
                                achance, amentor, anewb, ajobhunt, researchout))
                        researchmeans.append(pylab.mean(researchlist))
                        outfile.write(str(researchsum/sim_runs) + "\n")
                        researchses.append(
                            pylab.std(researchlist) / math.sqrt(sim_runs))

    datafile.close()
    gemfile.close()
    meansfile.close()
    outfile.close()

# ------------------------------------------------------------
# Run gem sweep
# ------------------------------------------------------------

def run_gem_sweep2(parms, seedrng):

    """
    Runs for sensitivity analysis using GEM-SA -- NO MENTORING
    """
    researchmeans = []
    researchses = []

    parms['mentored_pdrs'] = 0
    parms['write_output'] = False

    sim_runs = 60
    gemfile = open(parms['prefix'] + GEMSA_DATA + ".txt", 'w')
    meansfile = open(parms['prefix'] + GEMSA_MEANS + ".txt", 'w')
    outfile = open(parms['prefix'] + GEMSA_OUTS + ".txt", 'w')

    postdocchancelist = CHANCE_L
    newbtimelist = TIME2_L
    jobhunttimelist = TIME2_L

    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, True, True, True, False)
    root_dir = parms['prefix']

    for achance in postdocchancelist:
        for anewb in newbtimelist:
            for ajobhunt in jobhunttimelist:
                parms['jobhunt_time'] = ajobhunt
                parms['newb_time'] = anewb
                parms['postdoc_chance'] = achance
                print("Trying postdoc chance: ", achance)
                print("Trying newb time: ", anewb)
                print("Trying jobhunting time: ", ajobhunt)
                researchlist = []
                researchsum = 0.0
                meansfile.write(_tab_str3n(achance, anewb, ajobhunt))
                for arun in range(0, sim_runs):
                    print(arun)
                    researchout = _run_sim(datafile, parms, root_dir,
                                           seedrng, sim=True)[2]
                    researchlist.append(researchout)
                    researchsum += researchout
                    print(researchout)
                    gemfile.write(
                        _tab_str4n(achance, anewb, ajobhunt, researchout))
                    researchmeans.append(pylab.mean(researchlist))
                    outfile.write(str(researchsum/sim_runs) + "\n")
                    researchses.append(
                        pylab.std(researchlist) / math.sqrt(sim_runs))

    datafile.close()
    gemfile.close()
    meansfile.close()
    outfile.close()

# ------------------------------------------------------------
# Run gem sweep
# ------------------------------------------------------------

def run_gem_sweep3(parms, seedrng):

    """
    Runs for sensitivity analysis using GEM-SA -- INCLUDING MENTORING SWITCH
    """
    researchmeans = []
    researchses = []

    parms['write_output'] = False

    sim_runs = 2 #20
    gemfile = open(parms['prefix'] + GEMSA_DATA + ".txt", 'w')
    meansfile = open(parms['prefix'] + GEMSA_MEANS + ".txt", 'w')
    outfile = open(parms['prefix'] + GEMSA_OUTS + ".txt", 'w')

    postdocchancelist = CHANCE_L
    mentoringbonuslist = BONUS2_L
    newbtimelist = TIME1_L
    jobhunttimelist = TIME1_L

    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, True, True, True, False)
    root_dir = parms['prefix']

    for achance in postdocchancelist:
        for amentor in mentoringbonuslist:
            for anewb in newbtimelist:
                for ajobhunt in jobhunttimelist:
                    parms['postdoc_chance'] = achance
                    parms['mentoring_bonus'] = amentor
                    parms['newb_time'] = anewb
                    parms['jobhunt_time'] = ajobhunt
                    print("Trying postdoc chance: ", achance)
                    print("Trying mentoring bonus: ", amentor)
                    print("Trying newb time: ", anewb)
                    print("Trying jobhunting time: ", ajobhunt)
                    researchlist = []
                    researchsum = 0.0
                    meansfile.write(_tab_str4n(achance, amentor, anewb, ajobhunt))
                    for arun in range(0, sim_runs):
                        print(arun,)
                        researchout = _run_sim(
                            datafile, parms, root_dir, seedrng, sim=True)[2]
                        researchlist.append(researchout)
                        researchsum += researchout
                        print(researchout)
                        gemfile.write(
                            _tab_str5n(
                                achance, amentor, anewb, ajobhunt, researchout))
                        researchmeans.append(pylab.mean(researchlist))
                        outfile.write(str(researchsum/sim_runs) + "\n")
                        researchses.append(
                            pylab.std(researchlist) / math.sqrt(sim_runs))

    datafile.close()
    gemfile.close()
    meansfile.close()
    outfile.close()

# ------------------------------------------------------------
# Run gem sweep
# ------------------------------------------------------------

def run_gem_sweep4(parms, seedrng):

    """
    Runs for sensitivity analysis using GEM-SA -- INCLUDING MENTORING SWITCH
    """
    researchmeans = []
    researchses = []

    parms['write_output'] = False

    sim_runs = 2 #20
    gemfile = open(parms['prefix'] + GEMSA_DATA + ".txt", 'w')
    meansfile = open(parms['prefix'] + GEMSA_MEANS + ".txt", 'w')
    outfile = open(parms['prefix'] + GEMSA_OUTS + ".txt", 'w')

    postdocchancelist = CHANCE_L
    mentoringbonuslist = BONUS2_L
    newbtimelist = TIME1_L
    jobhunttimelist = TIME1_L

    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, True, True, True, False)
    root_dir = parms['prefix']

    for achance in postdocchancelist:
        for amentor in mentoringbonuslist:
            for anewb in newbtimelist:
                for ajobhunt in jobhunttimelist:
                    parms['postdoc_chance'] = achance
                    parms['mentoring_bonus'] = amentor
                    parms['newb_time'] = anewb
                    parms['jobhunt_time'] = ajobhunt
                    print("Trying postdoc chance: ", achance)
                    print("Trying mentoring bonus: ", amentor)
                    print("Trying newb time: ", anewb)
                    print("Trying jobhunting time: ", ajobhunt)
                    redundantlist = []
                    redundantsum = 0.0
                    meansfile.write(_tab_str4n(achance, amentor, anewb, ajobhunt))
                    for arun in range(0, sim_runs):
                        print(arun,)
                        #sim = _run_sim(datafile, parms, root_dir, seedrng, sim=True)
                        redundout = _run_sim(
                            datafile, parms, root_dir, seedrng, sim=True)[3]
                        redundantlist.append(redundout)
                        redundantsum += redundout
                        print(redundout)
                        gemfile.write(
                            _tab_str5n(
                                achance, amentor, anewb, ajobhunt, redundout))
                        researchmeans.append(pylab.mean(redundantlist))
                        outfile.write(str(redundantsum/sim_runs) + "\n")
                        researchses.append(
                            pylab.std(redundantlist) / math.sqrt(sim_runs))

    datafile.close()
    gemfile.close()
    meansfile.close()
    outfile.close()

# ------------------------------------------------------------
# Run gem sweep
# ------------------------------------------------------------

def run_gem_sweep5(parms, seedrng):

    """
    Runs for sensitivity analysis using GEM-SA -- INCLUDING MENTORING SWITCH
    """
    researchmeans = []
    researchses = []

    parms['write_output'] = False

    sim_runs = 1 #20
    gemfile = open(parms['prefix'] + GEMSA_DATA + ".txt", 'w')
    meansfile = open(parms['prefix'] + GEMSA_MEANS + ".txt", 'w')
    outfile = open(parms['prefix'] + GEMSA_OUTS + ".txt", 'w')

    mentoringbonuslist = BONUS1_L
    postdocchancelist = CHANCE_L
    startinggrantlist = GRANT_L
    yearlyincreaselist = INCREASE_L
    #mentoringbonuslist = BONUS2_L
    #newbtimelist = TIME1_L
    #jobhunttimelist = TIME1_L

    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, True, True, True, False)
    root_dir = parms['prefix']

    for achance in postdocchancelist:
        for amentor in mentoringbonuslist:
            for astartgrant in startinggrantlist:
                for anincrease in yearlyincreaselist:
                    parms['postdoc_chance'] = achance
                    parms['mentoring_bonus'] = amentor
                    parms['starting_grant_fund'] = astartgrant
                    parms['yearly_increase'] = anincrease
                    redundantlist = []
                    redundantsum = 0.0
                    meansfile.write(
                        _tab_str4n(achance, amentor, astartgrant, anincrease))
                    for arun in range(0, sim_runs):
                        print(arun,)
                        #sim = _run_sim(
                        #parms, root_dir, seedrng, sim=True)
                        redundout = _run_sim(
                            datafile, parms, root_dir, seedrng, sim=True)[3]
                        redundantlist.append(redundout)
                        redundantsum += redundout
                        print(redundout)
                        gemfile.write(
                            _tab_str5n(achance, amentor, astartgrant,
                                       anincrease, redundout))
                        researchmeans.append(pylab.mean(redundantlist))
                        outfile.write(str(redundantsum/sim_runs) + "\n")
                        researchses.append(
                            pylab.std(redundantlist) / math.sqrt(sim_runs))

    datafile.close()
    gemfile.close()
    meansfile.close()
    outfile.close()

# ------------------------------------------------------------
# Run gem sweep
# ------------------------------------------------------------

def run_gem_sweep6(parms, seedrng):

    """
    Runs for sensitivity analysis using GEM-SA --
    INCLUDING limited funding switch
    """
    researchmeans = []
    researchses = []

    parms['write_output'] = False

    sim_runs = 1 #20
    gemfile = open(parms['prefix'] + GEMSA_DATA + ".txt", 'w')
    meansfile = open(parms['prefix'] + GEMSA_MEANS + ".txt", 'w')
    outfile = open(parms['prefix'] + GEMSA_OUTS + ".txt", 'w')

    #startinggrantlist = GRANT_L
    #yearlyincreaselist = INCREASE_L
    #mentoringbonuslist = BONUS2_L
    #newbtimelist = TIME1_L
    jobhunttimelist = TIME2_L
    limitfunding = [True, False]
    mentoringbonuslist = BONUS1_L
    postdocchancelist = CHANCE_L

    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, True, True, True, False)
    root_dir = parms['prefix']

    for achance in postdocchancelist:
        for amentor in mentoringbonuslist:
            for ajobhunt in jobhunttimelist:
                for afund in limitfunding:
                    parms['postdoc_chance'] = achance
                    parms['mentoring_bonus'] = amentor
                    parms['jobhunt_time'] = ajobhunt
                    parms['limited_funding'] = afund
                    redundantlist = []
                    redundantsum = 0.0
                    meansfile.write(_tab_str4n(achance, amentor, ajobhunt, afund))
                    for arun in range(0, sim_runs):
                        print(arun,)
                        # sim = _run_sim(datafile, parms, root_dir, seedrng,
                        # sim=True)
                        redundout = _run_sim(
                            datafile, parms, root_dir, seedrng, sim=True)[3]
                        redundantlist.append(redundout)
                        redundantsum += redundout
                        print(redundout)
                        gemfile.write(
                            _tab_str5n(
                                achance, amentor, ajobhunt, afund, redundout))
                        researchmeans.append(pylab.mean(redundantlist))
                        outfile.write(str(redundantsum/sim_runs) + "\n")
                        researchses.append(
                            pylab.std(redundantlist) / math.sqrt(sim_runs))

    datafile.close()
    gemfile.close()
    meansfile.close()
    outfile.close()

# ------------------------------------------------------------
# Run a single run
# ------------------------------------------------------------

def run_once(parms, seedrng):

    """
    Run a single simulation
    """
    dfile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(dfile, False, False, False, False)
    _run_sim(dfile, parms, parms['prefix'], seedrng, sim=True)
    dfile.close()

# ------------------------------------------------------------
# Run partial *** needs to be looked at ***
# ------------------------------------------------------------

# Disable groupby checking in Pylint
# pylint: disable=E1101

def run_partial(rootdir):

    """
    Run partial -- needs to be looked at for functionality
    """
    plt.clf()
    fields = PDRGROW_L + RQPDRS_L + TOTROSACK_L
    dataf = pd.read_csv(rootdir + ROI_TEST + '.csv', usecols=fields)
    deviation = dataf.groupby(PDRGROW_L + RQPDRS_L).std()
    means = dataf.groupby(PDRGROW_L + RQPDRS_L).mean()
    print(means)
    print(deviation)
    scenariolist = ['Growing Pop', 'NoRQnoM', 'noRQM', 'RQnoM', 'RQM']

    # Total RO plot
    fig = _plot_fig(
        False, means[TOTRO], 'Output in Postdoc Scenarios',
        deviation.Total_RO, PALETTE1_L, "Scenarios", "Research Output",
        scenariolist, rootdir, MEMRUN_TOT, False)

    # Mean Output plot
    _plot_fig(
        fig, means[MEANROALL], 'Mean Output in Postdoc Scenarios',
        deviation.Mean_RO_all, PALETTE1_L, "Scenarios", "Mean Research Output",
        scenariolist, rootdir, MEMRUN_MEAN)

    # ROI in Postdoc Scenarios plot
    _plot_fig(
        fig, means['ROI_noPDRs'], 'ROI in Postdoc Scenarios',
        deviation.ROI_no_PDRs, 'g', "Scenarios", "Return on Investment",
        scenariolist, rootdir, MEMRUN_ROI)

    datafr = dataf[dataf.growing_pop != 1]
    deviationr = datafr.groupby([USEPDRS_L, RQPDRS_L]).std()
    meansr = datafr.groupby([USEPDRS_L, RQPDRS_L]).mean()
    scenariolistr = ['NoRQnoM', 'noRQM', 'RQnoM', 'RQM']

    # Redundancies in Postdoc Scenarios plot
    _plot_fig(
        fig, meansr[TOTSACK], 'Redundancies in Postdoc Scenarios',
        deviationr.Total_Sacked, 'r', "Scenarios", "Redundancies", scenariolistr,
        rootdir, MEMRUN_SACK)


# ------------------------------------------------------------
# Run sweep
# ------------------------------------------------------------

def run_sweep1(parms, seedrng):

    """
    Description of sweep
    """
    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, False, True, True, True)
    root_dir = parms['prefix']
    run_alife_sweep1(datafile, parms, parms['prefix'], seedrng)
    datafile.close()
    plt.clf()
    fields = FUNDLIM_L + TOTROSACK_L
    fundinglist = [True, False]

    # Fill data frame from CSV file
    dataf = pd.read_csv(root_dir + ROI_TEST + '.csv', usecols=fields)
    deviation = dataf.groupby(FUNDLIM).std()
    means = dataf.groupby(FUNDLIM).mean()
    print(means)
    print(deviation)

    # Limited Funding vs Output plot
    fig = _plot_fig(
        False, means[TOTRO], 'Limited Funding vs Output',
        deviation.Total_RO, 'c', "Funding Limit", "Research Output",
        fundinglist, root_dir, PRORUN_TOT, False)

    # Limited Funding vs Mean Output plot
    _plot_fig(
        fig, means['Mean_R0_all'], 'Limited Funding vs Mean Output',
        deviation.Mean_RO_all, 'm', "Funding Limit", "Mean Research Output",
        fundinglist, root_dir, PRORUN_MEAN)

    # Limited Funding vs ROI plot
    _plot_fig(
        fig, means[ROINOPDRS], 'Limited Funding vs ROI',
        deviation.ROI_no_PDRs, 'g', "Funding Limit", "Return on Investment",
        fundinglist, root_dir, PRORUN_ROI)

    plt.clf()
    #datafr = dataf[dataf.promo_chance != 1]
    #meansr = datafr.groupby(PROCHAN_L).mean()
    #deviationr = datafr.groupby(PROCHAN_L).std()

    # Limited Funding vs Redundancies plot
    _plot_fig(
        fig, means[TOTSACK], 'Limited Funding vs Redundancies',
        deviation.Total_Sacked, 'r', "Funding Limit", "Redundancies",
        fundinglist, root_dir, PRORUN_SACK, False)


# ------------------------------------------------------------
# Run sweep
# ------------------------------------------------------------

def run_sweep2(parms, seedrng):

    """
    Description of sweep
    """
    # Alife XV promo chance sweep
    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, False, True, True, False)
    root_dir = parms['prefix']
    run_alife_sweep2(datafile, parms, parms['prefix'], seedrng)
    datafile.close()
    plt.clf()

    fields = PROCHAN_L + TOTROSACK_L
    dataf = pd.read_csv(root_dir + ROI_TEST + '.csv', usecols=fields)
    deviation = dataf.groupby(PROCHAN).std()
    means = dataf.groupby(PROCHAN).mean()
    promochancelist = [15, 25, 50, 75, 100]
    print(means)
    print(deviation)

    # Promotion chance vs Output plot
    fig = _plot_fig(
        False, means[TOTRO], 'Promotion Chance vs Output',
        deviation.Total_RO, 'c', "Promotion Chance", "Research Output",
        promochancelist, root_dir, PRORUN_TOT)

    # Promotion chance vs Mean Output plot
    _plot_fig(
        fig, means[MEANROALL], 'Promotion Chance vs Mean Output',
        deviation.Mean_RO_all, 'm', "Promotion Chance", "Mean Research Output",
        promochancelist, root_dir, PRORUN_MEAN)

    # Promotion chance vs Mean Output plot
    _plot_fig(
        fig, means[ROINOPDRS], 'Promotion Chance vs ROI',
        deviation.ROI_no_PDRs, 'g', "Promotion Chance", "Return on Investment",
        promochancelist, root_dir, PRORUN_ROI)

    datafr = dataf[dataf.promo_chance != 1]
    deviationr = datafr.groupby(PROCHAN_L).std()
    meansr = datafr.groupby(PROCHAN_L).mean()

    # Redundancies vs Promotion chance plot
    _plot_fig(
        fig, meansr[TOTSACK], 'Redundancies vs Promotion Chance',
        deviationr.Total_Sacked, 'r', "Promotion Chance", "Redundancies",
        promochancelist, root_dir, PRORUN_SACK)


# ------------------------------------------------------------
# Run sweep
# ------------------------------------------------------------

def run_sweep3(parms, seedrng):

    """
    Description of sweep
    """
    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, False, True, True, True)
    root_dir = parms['prefix']
    run_alife_sweep3(datafile, parms, parms['prefix'], seedrng)
    datafile.close()
    plt.clf()
    fields = FUNDINC_L + TOTROSACK_L
    dataf = pd.read_csv(root_dir + ROI_TEST + '.csv', usecols=fields)
    deviation = dataf.groupby(FUNDINC).std()
    fundinglist = [2, 2.5, 3, 3.5, 4, 4.5, 5]
    means = dataf.groupby(FUNDINC).mean()
    print(means)
    print(deviation)

    # Funding Increases vs Output plot
    fig = _plot_fig(
        False, means[TOTRO], 'Funding Increases vs Output',
        deviation.Total_RO, 'c', "Funding Increase", "Research Output",
        fundinglist, root_dir, INCRUN_TOT, False)

    # Funding Increases vs Mean Output
    _plot_fig(
        fig, means[MEANROALL], 'Funding Increases vs Mean Output',
        deviation.Mean_RO_all, 'm', "Funding Increase", "Mean Research Output",
        fundinglist, root_dir, INCRUN_MEAN)

    # Funding Increases vs ROI
    _plot_fig(
        fig, means[ROINOPDRS], 'Funding Increases vs ROI',
        deviation.ROI_no_PDRs, 'g', "Funding Increase", "Return on Investment",
        fundinglist, root_dir, INCRUN_ROI)

    #datafr = dataf[dataf.promo_chance != 1]
    #meansr = datafr.groupby(PROCHAN_L).mean()
    #deviationr = datafr.groupby(PROCHAN_L).std()

    # Funding Increases vs Redundancies
    _plot_fig(
        fig, means[TOTSACK], 'Funding Increases vs Redundancies',
        deviation.Total_Sacked, 'r', "Funding Increase", "Redundancies",
        fundinglist, root_dir, INCRUN_SACK)

# ------------------------------------------------------------
# Run sweep
# ------------------------------------------------------------

def run_sweep4(parms, seedrng):

    """
    Description of sweep
    """
    # Parameter sweep, baseline vs limited vs unlimited funding
    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, False, True, True, True)
    root_dir = parms['prefix']
    run_alife_sweep4(datafile, parms, parms['prefix'], seedrng)
    datafile.close()
    plt.clf()
    fields = PDRGROW_L + TOTROSACK_L
    dataf = pd.read_csv(root_dir + ROI_TEST + '.csv', usecols=fields)
    deviation = dataf.groupby(PDRGROW_L).std()
    means = dataf.groupby(PDRGROW_L).mean()
    print(means)
    print(deviation)
    scenariolist = ['Growing Pop', 'Limited funding']

    # Output Comparison: Baseline vs Limited Funding
    fig = _plot_fig(
        False, means[TOTRO], 'Output Comparison: Baseline vs Limited Funding',
        deviation.Total_RO, PALETTE2_L, "Scenarios", "Research Output",
        scenariolist, root_dir, MEMRUN_TOT, False)

    # Mean Output Comparison: Baseline vs Limited Funding
    _plot_fig(
        fig, means[TOTSACK],
        'Mean Output Comparison: Baseline vs Limited Funding',
        deviation.Mean_RO_all, PALETTE3_L, "Scenarios", "Main Research Output",
        scenariolist, root_dir, MEMRUN_MEAN)

    # ROI in Limited-Funding Postdoc Scenarios
    _plot_fig(
        fig, means[ROINOPDRS], 'ROI in Limited-Funding Postdoc Scenarios',
        deviation.ROI_no_PDRs, 'g', "Scenarios", "Return on Investment",
        scenariolist, root_dir, MEMRUN_ROI)

    # Redundancies: Baseline vs Limited Funding
    _plot_fig(
        fig, means[TOTSACK], 'Redundancies: Baseline vs Limited Funding',
        deviation.ROI_no_PDRs, 'r', "Scenarios", "Redundancies",
        scenariolist, root_dir, INCRUN_SACK)

# ------------------------------------------------------------
# Run sweep
# ------------------------------------------------------------

def run_sweep5(parms, seedrng):

    """
    Description of sweep
    """
    # Parameter sweep, baseline vs limited vs unlimited funding
    datafile = open(parms['prefix'] + ROI_TEST + '.csv', 'w')
    _headers(datafile, False, True, True, True)
    root_dir = parms['prefix']
    run_alife_sweep5(datafile, parms, parms['prefix'], seedrng)
    datafile.close()
    plt.clf()
    fields = PDRGROW_L + TOTROSACK_L + FUNDLIM_L
    dataf = pd.read_csv(root_dir + ROI_TEST + '.csv', usecols=fields)
    deviation = dataf.groupby(PDRGROW_L + FUNDLIM_L).std()
    means = dataf.groupby(PDRGROW_L + FUNDLIM_L).mean()
    print(means)
    print(deviation)
    scenariolist = ['GP-LF', 'PDR-LF', 'GP-UF', 'PDR-UF']

    # Output Comparison: Baseline vs Limited Funding
    fig = _plot_fig(
        False, means[TOTRO], 'Output Comparison: Baseline vs Limited Funding',
        deviation.Total_RO, PALETTE4_L, "Scenarios", "Research Output",
        scenariolist, root_dir, MEMRUN_TOT, False)

    # Mean Output Comparison: Baseline and Funding Scenarios
    _plot_fig(
        fig, means[TOTSACK],
        'Mean Output Comparison: Baseline and Funding Scenarios',
        deviation.Mean_RO_all, PALETTE4_L, "Scenarios", "Mean Researhc Output",
        scenariolist, root_dir, MEMRUN_MEAN)

    # Mean Output Comparison: Baseline and Funding Scenarios
    _plot_fig(
        fig, means[ROINOPDRS], 'ROI by Scenario',
        deviation.ROI_nno_PDRs, 'g', "Scenarios", "Mean Research Output",
        scenariolist, root_dir, MEMRUN_MEAN)

    # ROI by Scenario
    _plot_fig(
        fig, means[ROINOPDRS], 'ROI by Scenario',
        deviation.ROI_nno_PDRs, 'g', "Scenarios", "Return on Investment",
        scenariolist, root_dir, MEMRUN_ROI)

    # Redundancies: Baseline and Funding Scenarios
    _plot_fig(
        fig, means[TOTSACK], 'Redundancies: Baseline and Funding Scenarios',
        deviation.Total_Sacked, 'r', "Scenarios", "Redundancies",
        scenariolist, root_dir, MEMRUN_SACK)

# pylint: enable=E1101

# ------------------------------------------------------------
# WCSS runs
# ------------------------------------------------------------

def run_wcss_sims(dfile, parms, base_prefix, seedrng):

    """
    Run basic simulations associated with the WCSS paper.
    """
    parms['runs'] = 10

    # run THERMOSTAT model
    parms['learning_type'] = 'thermostat'
    parms['prefix'] = base_prefix + "thermostat/"
    print('thermostat')
    _run_sim(dfile, parms, parms['prefix'], seedrng)

    # run MEMORY A model ("bad" parameters)
    parms['learning_type'] = 'memory'
    parms['prefix'] = base_prefix + "memory_A/"
    parms['prob_reentry'] = 0.05
    parms['run_length'] = 5
    print('memory A')
    _run_sim(dfile, parms, parms['prefix'], seedrng)

    # run MEMORY B model ("good" parameters)
    parms['learning_type'] = 'memory'
    parms['prefix'] = base_prefix + "memory_B/"
    parms['prob_reentry'] = 0.02
    parms['run_length'] = 3
    print('memory B')
    _run_sim(dfile, parms, parms['prefix'], seedrng)

    # run FIXED model
    parms['prefix'] = base_prefix + "fixed/"
    print('fixed')
    _run_sim(dfile, parms, parms['prefix'], seedrng, sim=False)

# ============================================================
# Public functions
# ============================================================

def init_params():

    """
    *** Description ***
    """
    parms = {}

    # simulation parameters
    parms['career_end'] = 60        # semesters before agent has to retire
    parms['fixed_time'] = 0.1       # ie, for legislated time alloc
    parms['growing_pop'] = 0        # use growing population, no postdocs
    parms['init_time'] = 0.5        # upper bound on init'l time_grant values
    parms['iterations'] = 100       # num of iterations to simulate
    parms['jobhunt_time'] = 0.3     # time spent job-hunting in final 2
                                    # semesters of contract
    parms['mentored_pdrs'] = 1      # do postdocs gain RQ due to mentoring
    parms['mentoring_bonus'] = 0.20 # bonus to RQ for promoted postdocs
                                    #(from being mentored/maturing)
    parms['newb_time'] = 0.4        # time spent being new postdoc/academic
    parms['pdr_rq_counts'] = 1      # does postdoc RQ count in promotions
    parms['pop_size'] = 100         # initial num of academic agents
    parms['postdoc_chance'] = 0.15  # chance for PDR promotion
    parms['prefix'] = DIR_PREFIX    # where to write output
    parms['random_seed'] = True     # whether to use random seed (or fixed)
    parms['runs'] = 1               # num of runs per param combination
    parms['seed'] = 1234            # seed to use (if random_seed==False)
    parms['use_postdocs'] = 1       # whether to include postdocs in sim
    parms['use_retirement'] = True
    parms['write_output'] = True    # whether or not to write output on runs

    # grant parameters
    parms['grant_bonus'] = 1.5      # G: bonus to res output from grants
    parms['grant_noise'] = 0.1      # std dev gaussian noise on grant quality
    parms['grant_pools'] = 1        # number of pools for grant evaluation
    parms['grant_proportion'] = 0.3 # P: prop of pop'n who can obtain grants
    parms['grant_slope'] = 2.0      # slope constant in tanh function
    parms['limited_funding'] = True # are grants limited?
    parms['manager_penalty'] = 0.00 # mgt time deducted from successful
                                    #grant applicants
    parms['research_slope'] = 2.0   # slope constant in tanh function
    parms['rq_counts'] = True       # is res_quality used in grant_quality?
    parms['starting_grant_fund'] = 30 # starting funds available,
                                    # 1 unit equals 1 funded project
    parms['weight_grant'] = 1.0     # weighting on grant qual rt track record
    parms['yearly_increase'] = 0.02 # % increase per step in ltd funding case

    # self learning parameters
    parms['learning_type'] = 'memory' # options 'thermostat', 'memory'
    parms['memory_size'] = 12       # num of memory steps to store
                                    # (NB: this is NOT window length!)
    parms['prob_reentry'] = 0.02    # E: probability of re-entering population
                                    # after dropping out
    parms['reentry_range'] = 0.2    # upper bound on re-entry time_grant
    parms['run_length'] = 12        # W: num of memory steps to consider
                                    # (THIS is window length!)
    parms['self_update_width'] = 0.1  # "learning rate"
    parms['self_update_width_fixed'] = True # fixed width learning steps?

    return parms

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
