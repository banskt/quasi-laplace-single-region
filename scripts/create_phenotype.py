import argparse
import random
import collections
import os
import numpy as np
import scipy.stats

SNPINFO_FIELDS = ['rsid', 'bp_location', 'alt_allele', 'ref_allele', 'studies']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

def parse_args():

    parser = argparse.ArgumentParser(description='Simulate phenotype given genotype')

    parser.add_argument('-dl', '--locidir',
                        type=str,
                        dest='locidir',
                        metavar='DIR',
                        help='path of the loci dosages directory')

    parser.add_argument('-ds', '--simdir',
                        type=str,
                        dest='simdir',
                        metavar='DIR',
                        help='path of the simulation directory')

    parser.add_argument('-fl', '--locusnames',
                        type=str,
                        dest='locusnames',
                        metavar='FILE',
                        help='file with prefix of locusnames to use')

    parser.add_argument('-st', '--studies',
                        nargs='*',
                        type=str,
                        dest='studies',
                        metavar='STR',
                        help='list of study names')

    parser.add_argument('-stn', '--samples',
                        nargs='*',
                        type=int,
                        dest='samples',
                        metavar='INT',
                        help='list of numbers of samples in studies')

    parser.add_argument('-sp', '--sampleprefix',
                        default='_QC_age',
                        type=str,
                        dest='sampleprefix',
                        metavar='STR',
                        help='the prefix of the input sample file named as {STUDY}{SAMPLEPREFIX}.sample')

    parser.add_argument('-o', '--outfile',
                        default='causal.snplist',
                        type=str,
                        dest='outfile',
                        metavar='FILE',
                        help='name of file for causal SNPs list')

    parser.add_argument('-p', '--prop',
                        default=5,
                        type=float,
                        dest='proportion',
                        metavar='FLOAT',
                        help='fraction of causal SNPs')

    parser.add_argument('-hg', '--heritability',
                        default=0.25,
                        type=float,
                        dest='sigma_herited_sq',
                        metavar='real',
                        help='genotypic variance (sigma^2) / narrow-sense heritability ')

    parser.add_argument('-k', '--prevalence',
                        default=0.3,
                        type=float,
                        dest='prevalence',
                        metavar='real',
                        help='disease prevalence')

    parser.add_argument('-t', '--simtype',
                        default='fixed',
                        type=str,
                        dest='simtype',
                        metavar='STR',
                        help='type of distribution to be used for effectsizes')

    parser.add_argument('--agesex', 
                        dest='use_agesex',
                        action='store_true',
                        help='whether to use age and sex as covariates')


    parser.add_argument('-cr', '--ccratio',
                        default=1.0,
                        type=float,
                        dest='ccratio',
                        metavar='real',
                        help='ratio of cases to controls, cannot be greater than 1 ')

    opts = parser.parse_args()
    return opts


def read_locus_prefixes(filename):
    locusprefixes = list()
    with open(filename, 'r') as mfile:
        for line in mfile:
            mline = line.split()
            locusprefixes.append(mline[0].strip())
    return locusprefixes

def read_snps(locidir, studies, locusprefixes):
    nstudy = len(studies)
    snpinfo = list()
    for study in studies:
        snpstudy = list()
        for locusprefix in locusprefixes:
            mapfile = os.path.join(locidir, study, locusprefix + '.map')
            snplocus = list()
            with open(mapfile, 'r') as mfile:
                for line in mfile:
                    mline = line.split()
                    this_snp = SnpInfo(rsid = mline[1],
                                       bp_location = int(mline[3]),
                                       alt_allele = mline[4],
                                       ref_allele = mline[5],
                                       studies = [])
                    snplocus.append(this_snp)
            snpstudy.append(snplocus)
        snpinfo.append(snpstudy)
    return snpinfo


def norm_binom(gt, f):
    gt = (gt - (2 * f)) / np.sqrt(2 * f * (1 - f))
    return gt


def read_genotype(gendir, snpinfo, nsample):
    ncol = nsample * 3 + 5
    nloci = len(snpinfo)
    genotypelist = list()
    for l, locusprefix in enumerate(locusprefixes):
        filename = os.path.join(gendir, locusprefix + '.gen')
        dosage = np.loadtxt(filename, usecols=range(5, ncol)).T
        nsnps = len(snpinfo[l])
        genotype = np.empty([nsample, nsnps], dtype=np.float64)
        for j in range(nsample):
            maj_ind = j * 3            # Allele AA, where A is the alternate allele
            het_ind = maj_ind + 1      # Allele AB, where B is the reference allele
            min_ind = het_ind + 1      # Allele BB
            genotype[j,:] = 2 * dosage[min_ind, :] + dosage[het_ind, :] # [AA, AB, BB] := [0, 1, 2]
        freq  = np.sum(genotype, axis=0) / (2 * genotype.shape[0])
        genotype = norm_binom(genotype, freq)
        genotypelist.append(genotype)
    return genotypelist

    
def get_snppool(snpinfo):
    nstudy = len(snpinfo)
    snppool = list()
    nloci = len(snpinfo[0]) # Assumes same number of loci in each study
    for locus in range(nloci):
        snplocus = list()
        for i in range(nstudy):
            snplocus += [snp for snp in snpinfo[i][locus] if snp not in snplocus]
        snppool.append(snplocus)
    return snppool


def select_causal_snps(studies, snpinfo, snppool, locus, prob):
    select = list()
    none_selected = True
    while none_selected:
        for i, snp in enumerate(snppool):
            mrand = random.uniform(0, 1)
            if prob > mrand:
                select.append(snp)
        if len(select) > 0:
            none_selected = False
        #else:
        #    select = list()
    for snp in select:
        present_in = [study for i, study in enumerate(studies) if snp in snpinfo[i][locus]]
        select[select.index(snp)] = snp._replace(studies = present_in)
    return select

# ================= Phenotype simulation ===================

opts = parse_args()
studies = opts.studies
samples = opts.samples
cprob = opts.proportion
insampleprefix = opts.sampleprefix
locidir = os.path.realpath(opts.locidir)
simdir  = os.path.realpath(opts.simdir)
sampledir = os.path.join(simdir, 'samples')
if not os.path.exists(sampledir):
    os.makedirs(sampledir)
effectfilename = 'snps.effectsize'
samplefilename = 'phenotypes.sample'
ccratio = opts.ccratio

locusprefixes = read_locus_prefixes(opts.locusnames)
nloci = len(locusprefixes)

# Read all SNP info
snpinfo = read_snps(locidir, studies, locusprefixes)
snppool = get_snppool(snpinfo)
print ("Read all SNP info")

# Select causal SNPs per loci and write in causal.snplist
causality = np.zeros(nloci)
causal_snps = list()
ncausal = 0
fname = os.path.join(sampledir, opts.outfile)
with open(fname, 'w') as mfile:
    for l, locus in enumerate(locusprefixes):
        snps = select_causal_snps(studies, snpinfo, snppool[l], l, cprob)
        ncausal += len(snps)
        if len(snps) > 0:
            causality[l] = 1
        causal_snps.append(snps)
        mfile.write('#--------------------\n')
        mfile.write('Causality: {:s} {:g}\n'.format(locus, causality[l]))
        for snp in causal_snps[l]:
            mfile.write('{:s} {:d}\t{:d} [{:s}]\n'.format(snp.rsid, snp.bp_location, len(snp.studies), ' '.join(snp.studies)))

print ("{:g} causal loci, {:d} causal SNPs, zmax = {:d}".format(np.sum(causality), ncausal, max([len(x) for x in causal_snps if x is not None])) )

beta = np.zeros(ncausal)
if opts.simtype == 'fixed':
    beta = np.ones(ncausal)
elif opts.simtype == 'bimodal':
    mean1 = 0.5
    mean2 = -0.5
    mvar = 0.2
    for i in range(ncausal):
        mrandom = np.random.uniform()
        if mrandom <= 0.5:
            beta[i] = np.random.normal(mean1, mvar)
        else:
            beta[i] = np.random.normal(mean2, mvar)
elif opts.simtype == 'studentsT':
    beta = np.random.standard_t(1, size = ncausal)
elif opts.simtype == 'normal':
    beta = np.random.normal(0, 1, size = ncausal)
else:
    beta = np.random.rand(ncausal)

beta *= np.sqrt( opts.sigma_herited_sq / np.sum(np.square(beta)) )

# Simulate phenotype for each study
for i, study in enumerate(studies):
    print ("Creating phenotype for cohort %s" % study)

    # File paths
    outdir = os.path.join(sampledir, study)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    gendir = os.path.join(locidir, study)
    insample = os.path.join(gendir, '{:s}{:s}.sample'.format(study, insampleprefix))
    outsample = os.path.join(outdir, samplefilename)

    # Read genotype
    nsample = samples[i]
    target_cases = int(nsample / 2.1)
    gt = read_genotype(gendir, snpinfo[i], nsample)
    gt_tot = np.concatenate(gt, axis=1)
    allrsid = [snp.rsid for locus in snpinfo[i] for snp in locus]  # flatten out the snpinfo.rsid

    fname = os.path.join(outdir, effectfilename)
    mfile = open(fname, 'w')
    gtmask = np.zeros(gt_tot.shape[1], dtype=bool)
    betamask = np.zeros(beta.shape[0], dtype=bool)
    count = 0
    for locus in range(nloci):
        if causality[locus] == 1:
            for snp in causal_snps[locus]:
                if snp.rsid in allrsid:
                    k = allrsid.index(snp.rsid)
                    gtmask[k] = True
                    betamask[count] = True
                    mfile.write("%s \t %g\n" % (allrsid[k], beta[count]))
                else:
                    mfile.write("No SNPs found with rsid %s\n" % snp.rsid)
                count += 1
    mfile.close()

    gt_causal = gt_tot[:, gtmask]
    beta_sel = beta[betamask]

    # Get disease liabilities
    geno_eff = np.sum((beta_sel * gt_causal), axis=1)
    rand_var = np.sqrt(1 - opts.sigma_herited_sq)
    cov_eff = np.zeros(geno_eff.shape[0])
    if opts.use_agesex:
        age = np.loadtxt(insample, skiprows=2, usecols=6)
        sex = np.loadtxt(insample, skiprows=2, usecols=5)
        cov_eff[np.where(sex == 1)] = age[np.where(sex==1)] * 0.049 - 2.691
        cov_eff[np.where(sex == 2)] = age[np.where(sex==2)] * 0.055 - 3.197
        cov_var_male    = 0.049 * 0.049 * np.var(age[np.where(sex==1)])
        cov_var_female  = 0.055 * 0.055 * np.var(age[np.where(sex==2)])
        rand_var_male   = np.sqrt(max(0, 1 - opts.sigma_herited_sq - cov_var_male))
        rand_var_female = np.sqrt(max(0, 1 - opts.sigma_herited_sq - cov_var_female))
        rand_eff = np.zeros(geno_eff.shape[0])
        rand_eff[np.where(sex == 1)] = np.random.normal(0, rand_var_male,   np.where(sex == 1)[0].shape[0])
        rand_eff[np.where(sex == 2)] = np.random.normal(0, rand_var_female, np.where(sex == 2)[0].shape[0])
        print ("Variance of the random term. Males: {:g} Females: {:g}\n".format(np.square(rand_var_male), np.square(rand_var_female)))
    else:
        rand_eff = np.random.normal(0, rand_var, geno_eff.shape[0])

    liability = geno_eff + cov_eff + rand_eff
    # Select liabilities above the threshold of normal distribution truncating at K (disease prevalence)
    cutoff = scipy.stats.norm.ppf(1 - opts.prevalence)
    cases = np.where(liability >= cutoff)[0]
    pheno = np.zeros(nsample, dtype=int)
    pheno[cases] = 1
    keepsample = np.ones(nsample, dtype=bool)
    ncase = int(cases.shape[0])
    nctrl = int(nsample - ncase)
    if ccratio < 1.0:
        ncase = int(nctrl * ccratio)
        caseprob = 1 - scipy.stats.norm.pdf(liability[cases], 0, 1)
        normed_caseprob = caseprob / np.sum(caseprob)
        choosecase = np.random.choice(cases, size=ncase, replace=False, p=normed_caseprob)
        droppedcase = np.array([x for x in cases if x not in choosecase])
        keepsample[droppedcase] = False

    # Copy the sample file with new phenotype
    k = 0
    with open(insample, 'r') as samfile, open(outsample, 'w') as outfile:
        line1 = samfile.readline()
        line2 = samfile.readline()
        outfile.write(line1)
        outfile.write(line2)
        for line in samfile:
            if keepsample[k]:
                newpheno = line.strip()[:-1] + '{:g}\n'.format(pheno[k])
                outfile.write(newpheno)
            else:
                newpheno = line.strip()[:-1] + 'NA\n'
                outfile.write(newpheno)
            k += 1

    print ("Artificial phenotypes generated.")
    print ("No. of cases:", ncase)
    print ("No.of controls:", nctrl)
    print ()
