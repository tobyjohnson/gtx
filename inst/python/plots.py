

# on GSK systems requires /GWD/appbase/tools/bin/python2.7

from optparse import OptionParser
## deprecated, but argparse requires python >= 2.7 does not work on GSK servers
import logging
import gzip
import os
import sys
import string
import re
import math
import datetime

###
### parse command line arguments and check
###

parser = OptionParser(description = 'usage: %prog OPTIONS')
parser.add_option('-r', '--results', help = 'results file containing p-values', metavar = 'FILE',
                  type = 'string', dest = 'results')
parser.add_option('-o', '--output', help = 'output directory', metavar = 'DIRECTORY',
                  type = 'string', dest = 'output')

(options, args) = parser.parse_args()

logging.basicConfig(stream=sys.stdout, format='%(levelname)s:%(message)s', level=logging.DEBUG)
logging.info('genome wide plotting v0.1 Toby.x.Johnson@gsk.com')

if options.results is None:
    logging.error('No results file specified; use --results option')
    sys.exit(1)

if options.output is None:
    logging.error('No output directory, use --output option')
    sys.exit(1)

pos = dict() # dict of lists to store positions
pval = dict() # dict of lists to store -log10 P values
for chridx in range(22):
    pos[str(chridx + 1)] = []
    pval[str(chridx + 1)] = []    
#Causes sort to fail when drawing Manhatten plot - need to figure out
#pos['X'] = []
#pval['X'] = []

try:
    input_fh = gzip.open(options.results, 'r')
except OSError:
    logging.error('Could not read file [ ' + options.results + ' ]')
    sys.exit(2)

input_fh.readline() # header
for inputline in input_fh:
    inputdata = inputline.strip().split()
    pos[inputdata[0]].append(int(inputdata[1]))
    try:
        pval[inputdata[0]].append(-math.log10(float(inputdata[5])))
    except ValueError:
        pval[inputdata[0]].append(0.)

input_fh.close()

logging.info('Starting numpy and matplotlib modules')
import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import scipy.stats as stats

logging.info('Building numpy arrays...')
midpt = dict()
nudgelen = 50000000
nudge = nudgelen
pmax = 8.0
alltuple = list()
for chrom in sorted(pos.keys(), key=int):
    if len(pos[chrom]) > 0:
        logging.info(' ...for chr' + chrom)
        pos[chrom] = np.array(pos[chrom])
        pos[chrom] += nudge
        midpt[chrom] = (pos[chrom].min() + pos[chrom].max())/2
        nudge = pos[chrom].max() + nudgelen
        pval[chrom] = np.array(pval[chrom])
        pmax = max(pmax, pval[chrom].max())
        alltuple.append(pval[chrom])

logging.info(' ...for all chromosomes')
allpval = np.sort(np.concatenate(alltuple))

# compute median by simple lookup since array is already sorted
if len(allpval) % 2:
    pmedian = allpval[(len(allpval)-1)/2]
else:
    pmedian = 0.5*(allpval[len(allpval)/2-1]+allpval[len(allpval)/2])
# lambda is a reserved word in python
gclambda = stats.chi2.isf(pow(10, -pmedian), 1)/stats.chi2.isf(0.5, 1)
logging.info('Median P = ' + ("%.4f" % pow(10, -pmedian)) + ' lambda = ' + ("%.4f" % gclambda))

logging.info('Generating QQ plot [ ' + options.output + '/qq.png' + ' ]')
timebegin = datetime.datetime.now()

plt.figure(figsize = [6, 6])
plt.axis([0, math.ceil(pmax + 0.5), 0, math.ceil(pmax + 0.5)])
#plt.setp(plt.gca(), 'xticks', midpt.values())
#plt.setp(plt.gca(), 'xticklabels', midpt.keys())
plt.tick_params(bottom = 'on', left = 'on', top = 'off', right = 'off', direction='out')
plt.hold(1) # do we need this?
plt.plot([0, math.ceil(pmax + 0.5)], [0, math.ceil(pmax + 0.5)], 'r-')
plt.text(0.5*math.ceil(pmax + 0.5), 0.1, 'lambda = ' + ("%.4f" % gclambda), ha = 'center', va = 'baseline')
this = plt.plot(-np.log10((np.array(range(len(allpval))) + 0.5)/len(allpval)),
                allpval[::-1], '.') # stride -1 to reverse
plt.setp(this, 'color', [0, 0, 0], 'marker', '.', 'markeredgewidth', 0)
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.savefig(options.output + '/qq.png', dpi=200, bbox_inches = 'tight')

timeend = datetime.datetime.now()
logging.info('Generating QQ plot took ' + str((timeend-timebegin).seconds) + ' seconds')

logging.info('Generating Manhattan plot [ ' + options.output + '/manhattan.png' + ' ]')
timebegin = datetime.datetime.now()

plt.figure(figsize = [10, 6])
plt.axis([0, nudge, 0, math.ceil(pmax + 0.5)])
plt.setp(plt.gca(), 'xticks', midpt.values())
plt.setp(plt.gca(), 'xticklabels', midpt.keys())
plt.tick_params(bottom = 'on', left = 'on', top = 'off', right = 'off', direction='out')
plt.hold(1) # do we need this?
for chrom in sorted(pos.keys(), key=int):
    this = plt.plot(pos[chrom], pval[chrom], '.')
    # hard coded alternating darkcyan/lightgrey
    if int(chrom) % 2:
        plt.setp(this, 'color', [0, .5, .5], 'marker', '.', 'markeredgewidth', 0)
    else:
        plt.setp(this, 'color', [.75, .75, .75], 'marker', '.', 'markeredgewidth', 0)

plt.xlabel('Genomic position by chromosome')
plt.ylabel('Association -log10(p)')
plt.savefig(options.output + '/manhattan.png', dpi=200, bbox_inches = 'tight')

timeend = datetime.datetime.now()
logging.info('Generating Manhattan plot took ' + str((timeend-timebegin).seconds) + ' seconds')

logging.info('Finished plotting')
sys.exit(0)
