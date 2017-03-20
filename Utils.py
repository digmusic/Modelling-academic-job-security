"""
--------------------------------------------------------------------------------
FILE: Utils.py
DESC: Utility functions for the simulation
AUTH: thorsilver
VERS: 1.0 March 2017
REQS: Python 3.x (version 3.6 used)
--------------------------------------------------------------------------------
"""
from __future__ import print_function # ensure Python 3.x print form known
import sys
import os
import math
import pylab

def create_dir(dir_name):

    """
    *** Needs description ***
    """
    # create results directory if it doesn't already exist
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)


def transpose(matr):

    """
    Transpose the rectangular two-dimentional matrix matr.
    """
    return [[matr[y][x]
             for y in range(len(matr))] for x in range(len(matr[0]))]


def get_mean(data):

    """
    *** Needs description ***
    """
    if len(data) > 0:
        return sum(data)/len(data)
    else:
        return 0


def get_mean_and_sum(data):

    """
    Takes a 2D array of data and calculates mean and sum of each sub-list.
    """
    return [get_mean(x) for x in data], [sum(x) for x in data]


def create_colour_dict(keymax=0):

    """
    Create a colour dictionary with an entry for each state.
    """
    colour_dict = {}
    colour_dict['none'] = 'white' # nodes in no group
    index = 0
    keys = range(0, keymax)     # *** was xrange ***
    prev_colour = 0
    for g_id in [x for x in keys]:
        cur_colour = (prev_colour+383) % keymax
        colour_dict[g_id] = '#%02x%02x%02x' \
                % hsl2rgb(float(cur_colour)/(len(keys)), 1.0, 0.5)
        prev_colour = cur_colour
        index += 1
    return colour_dict


def get_red_green_colour(tomatch):

    """
    Get a colour in the red-green colour range corresponding to tomatch.

    tomatch in range [0.0, 1.0]
    0.0 -> red
    0.5 -> yellow
    1.0 -> green
    """
    return '#%02x%02x%02x' % hsl2rgb(tomatch/4.0, 1.0, 0.5)


def create_plot_colour_dict():

    """
    *** Needs description ***
    """
    colour_dict = {}
    colour_dict['none'] = 'white'
    colour_dict[0] = 'teal'
    colour_dict[1] = 'limegreen'
    colour_dict[2] = 'indigo'
    colour_dict[3] = 'orange'
    colour_dict[4] = 'plum'
    colour_dict[5] = 'flerp'
    colour_dict[6] = 'derp'
    return colour_dict


def hsl2rgb(hue, sat, light):

    """
    *** Needs description ***
    """
    if sat == 0.0:                # HSL values = [0-1]
        red = 255*light               # RGB results = [0-1]
        grn = 255*light
        blu = 255*light
    else:
        if light < 0.5:
            upper = light * (1 + sat)
        else:
            upper = (light + sat) - (light * sat)
        lower = 2 * light - upper

        red = rounded(255*hue_2_rgb(lower, upper, hue + (1.0/3.0)))
        grn = rounded(255*hue_2_rgb(lower, upper, hue))
        blu = rounded(255*hue_2_rgb(lower, upper, hue - (1.0/3.0)))
    return (red, grn, blu)


def hue_2_rgb(lower, upper, thresh):         # hue_2_rgb

    """
    *** Needs description ***
    """
    if thresh < 0:
        thresh += 1.0
    if thresh > 1:
        thresh -= 1.0
    if (6 * thresh) < 1:
        return lower + ((upper - lower) * 6 * thresh)
    if (2 * thresh) < 1:
        return upper
    if (3 * thresh) < 2:
        return lower + ((upper - lower) * 6 * ((2.0/3.0) - thresh))
    return lower


def rounded(floatval):

    """
    *** Needs description ***
    """
    return int(math.floor(floatval+0.5))


def write_plot(
        x_data, y_data, outfile, title_text, xlabel_text, ylabel_text, colours,
        labels=None, ylim=None, typ='lin', marker='', keys=None, linew=2):

    """
    A general x-y plotting function.
    """

    assert len(x_data) == len(y_data)
    pylab.clf()
    for indx in range(len(x_data)):
        if typ == 'lin':
            plot_string = 'pylab.plot('
        elif typ == 'log':
            plot_string = 'pylab.loglog('
        elif typ == 'semilogy':
            plot_string = 'pylab.semilogy('
        elif typ == 'semilogx':
            plot_string = 'pylab.semilogx('
        else:
            print("NetIO_plots::write_plot - unknown plot type.")
            sys.exit()

        if keys != None:
            colour = get_red_green_colour(keys[indx])
        else:
            colour = colours[indx]

        plot_string += "x_data["+str(indx)+"], y_data["+str(indx)+"], " +\
                "'" + colour + "'" + ", marker='"+marker+"'" +\
                ", linewidth=%d)" % (linew)

        #print(plot_string)

        exec(plot_string)

    pylab.xlabel(xlabel_text)
    pylab.ylabel(ylabel_text)
    pylab.title(title_text)
    if labels != None:
        lgd = pylab.legend(
            labels, loc=9, bbox_to_anchor=(0.5, -0.1), ncol=len(labels))
    if ylim != None:
        pylab.ylim(ylim[0], ylim[1])
    if labels != None:
        pylab.savefig(
            outfile+".pdf", format='pdf', additional_artists=lgd,
            bbox_inches="tight")
        pylab.savefig(
            outfile+".png", format='png', additional_artists=lgd,
            bbox_inches="tight")
    else:
        pylab.savefig(outfile+".pdf", format='pdf')
        pylab.savefig(outfile+".png", format='png')


def write_data(data, outfile, sep='\n'):

    """
    *** Needs description ***
    """
    out = open(outfile, 'w')
    for item in data:
        out.write(str(item) + sep)
    out.flush()
    out.close()


def write_data_2d(data, outfile, sep1=',', sep2='\n'):

    """
    *** Needs description ***
    """
    out = open(outfile, 'w')
    for outer in data:
        for inner in outer:
            out.write(str(inner) + sep1)
        out.write(sep2)
    out.flush()
    out.close()
