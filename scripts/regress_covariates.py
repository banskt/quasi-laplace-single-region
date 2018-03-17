import numpy as np
import os
import argparse
from sklearn.linear_model import LinearRegression

def parse_args():

    parser = argparse.ArgumentParser(description='Regress the phenotype with the covariates')

    parser.add_argument('-p', '--pheno',
                        type=str,
                        dest='phenofile',
                        metavar='FILE',
                        help='phenotype file in BIMBAM format')

    parser.add_argument('-c', '--cov',
                        type=str,
                        dest='covfile',
                        metavar='FILE',
                        help='covariates file in BIMBAM format with intercepts')

    parser.add_argument('-o', '--out',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='output file in BIMBAM format containing regressed phenotype')

    opts = parser.parse_args()
    return opts

opts = parse_args()
phenofile = opts.phenofile
covfile   = opts.covfile
outfile   = opts.outfile

pheno = np.loadtxt(phenofile).reshape(-1,1)
cov = np.loadtxt(covfile)
lm = LinearRegression(fit_intercept=False, normalize=False, copy_X=True, n_jobs=-1)
lm.fit(cov, pheno)
regressed = pheno - lm.predict(cov)
np.savetxt(outfile, regressed, fmt='%g')
