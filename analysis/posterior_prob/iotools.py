import numpy as np
import os
from get_plotvals import LDResult

def read_simres(simdir, key, locusprefixes, causal_rsids, zmax, muvar):
    res = None
    blore_path   = 'blore/meta_without_feature/zmax{:d}_mu{:s}_pi0.01_sig0.01/blore_meta_res/{:s}.gen.res'
    probit_path  = 'pimass/c{:d}_1e6_probit/output/{:s}.mcmc.txt'
    linear_path  = 'pimass/c{:d}_1e6_linear/output/{:s}.mcmc.txt'
    finemap_path = 'finemap/c{:d}/{:s}.snp'
    if key == 'blore':   res = blore  (simdir, blore_path,   locusprefixes, causal_rsids, zmax, muvar)
    if key == 'finemap': res = finemap(simdir, finemap_path, locusprefixes, causal_rsids, zmax)
    if key == 'probit':  res = pimass (simdir, probit_path,  locusprefixes, causal_rsids, zmax)
    if key == 'linear':  res = pimass (simdir, linear_path,  locusprefixes, causal_rsids, zmax)
    return res


def read_locusprefixes(locusfile):
    with open(locusfile, 'r') as mfile:
        locusprefixes = mfile.readlines()
    locusprefixes = [x.strip() for x in locusprefixes]
    return locusprefixes


def read_causal_rsids(causal_snps_file):
    causal_rsid = dict()
    with open(causal_snps_file, 'r') as mfile:
        rsidlist = list()
        locus_name = ''
        for line in mfile:
            if line.startswith('#'):
                causal_rsid[locus_name] = rsidlist
            else:
                mline =  line.split()
                if mline[0].startswith('Causality'):
                    locus_name = mline[1]
                    rsidlist = list()
                else:
                    rsidlist.append(mline[0])
        causal_rsid[locus_name] = rsidlist
    return causal_rsid


def blore(simdir, filepath, locusprefixes, causal_rsids, zmax, muvar):
    thisres = list()
    for locus in locusprefixes:
        locusres = list()
        outfile = os.path.join(simdir, filepath.format(zmax, muvar, locus))
        causals = causal_rsids[locus]
        with open(outfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mline_split = mline.split()
                rsid = mline_split[0]
                prob = float(mline_split[4])
                causality = 1 if rsid in causals else 0
                mres = LDResult(locus = locus,
                                rsid = rsid,
                                stat = prob,
                                ld = 1,
                                causality = causality)
                locusres.append(mres)
        thisres.append(locusres)
    return thisres


def pimass(simdir, filepath, locusprefixes, causal_rsids, zmax):
    thisres = list()
    for locus in locusprefixes:
        locusres = list()
        outfile = os.path.join( simdir, filepath.format(zmax, locus))
        causals = causal_rsids[locus]
        with open(outfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mline_split = mline.split()
                rsid = mline_split[0].strip()
                prob = float(mline_split[3])
                causality = 1 if rsid in causals else 0
                mres = LDResult(locus = locus,
                                rsid = rsid,
                                stat = prob,
                                ld = 1,
                                causality = causality)
                locusres.append(mres)
        thisres.append(locusres)
    return thisres


def finemap(simdir, filepath, locusprefixes, causal_rsids, zmax):
    thisres = list()
    for locus in locusprefixes:
        locusres = list()
        outfile = os.path.join( simdir, filepath.format(zmax, locus))
        causals = causal_rsids[locus]
        with open(outfile, 'r') as mfile:
            next(mfile)
            for mline in mfile:
                mlinesplit = mline.split()
                rsid = mlinesplit[1]
                prob = float(mlinesplit[3])
                causality = 1 if rsid in causals else 0
                mres = LDResult(locus = locus,
                                rsid = rsid,
                                stat = prob,
                                ld = 1,
                                causality = causality)
                locusres.append(mres)
        thisres.append(locusres)
    return thisres
