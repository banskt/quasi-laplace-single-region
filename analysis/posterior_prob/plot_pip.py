import numpy as np
import os
import argparse
import collections

import stylesheet
import precision_recall_scores as prs
import iotools
import get_plotvals

def parse_args():

    parser = argparse.ArgumentParser(description='Plot precision and recall of PIPs')

    parser.add_argument('--outfile',
                        type=str,
                        dest='outfile',
                        metavar='STR',
                        help='path of the output file')

    parser.add_argument('--whichplot',
                        nargs='*',
                        type=str,
                        dest='whichplot',
                        metavar='STR',
                        help='list of methods, allowed values: blore, finemap, probit, linear')

    parser.add_argument('--basedir',
                        type=str,
                        dest='basedir',
                        metavar='STR',
                        help='where to find the simulation directories')

    parser.add_argument('--locusfile',
                        type=str,
                        dest='locusfile',
                        metavar='STR',
                        help='file containing the locus prefixes to be analyzed')

    parser.add_argument('--start',
                        type=int,
                        dest='startsim',
                        metavar='INT',
                        help='starting simulation number')

    parser.add_argument('--end',
                        type=int,
                        dest='endsim',
                        metavar='INT',
                        help='last simulation number')

    parser.add_argument('--cmax',
                        type=int,
                        dest='cmax',
                        metavar='INT',
                        help='maximum number of causal SNPs allowed')


    parser.add_argument('--credible',
                        dest='use_credible',
                        action='store_true',
                        help='if set, ranks each locus separately')


    opts = parser.parse_args()
    return opts


opts = parse_args()
locusfile = os.path.abspath(opts.locusfile)
outfile = os.path.abspath(opts.outfile)
basedir = os.path.abspath(opts.basedir)

if not os.path.exists(os.path.dirname(outfile)):
    os.makedirs(os.path.dirname(outfile))

locusprefixes = iotools.read_locusprefixes(locusfile)

res = collections.defaultdict(lambda:0)
for key in opts.whichplot:
    res[key] = list()

for sim in range(opts.startsim, opts.endsim + 1):
    simname = 'sim{:03d}'.format(sim)
    print ('Reading {:s}'.format(simname))
    simdir = os.path.join( basedir, simname)
    causal_snps_file = os.path.join(simdir, 'samples', 'causal.snplist')
    causal_rsids = iotools.read_causal_rsids(causal_snps_file)
    for key in opts.whichplot:
        thisres = iotools.read_simres(simdir, key, locusprefixes, causal_rsids, opts.cmax, '0')
        res[key].append(thisres)

nmax = 0
plotvals = collections.defaultdict(lambda:0)
for key in opts.whichplot:
    if opts.use_credible:
        data = list()
        for x in res[key]:
            data += x
    else:
        data = [[y for z in x for y in z] for x in res[key]]
    maxlen = max([len(x) for x in data])
    if maxlen > nmax:
        nmax = maxlen
    plotvals[key] = get_plotvals.precision_recall_threshold(data)



xlim = [0, int(0.15 * nmax)]
ylim = [0, 0.9]
xticks = None
yticks = np.arange(0.1, 0.9, 0.1)

stylesheet.save_prcplot(outfile, plotvals, xlim, ylim, xticks, yticks)
