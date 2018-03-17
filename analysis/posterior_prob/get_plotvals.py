import numpy as np
import collections
import precision_recall_scores as roc

INFO_FIELDS = ['locus', 'rsid', 'stat', 'ld', 'causality']
class LDResult(collections.namedtuple('_LDResult', INFO_FIELDS)):
    __slots__ = ()

''' Calculate the average recall and average precision at each threshold.
    The average is calculated by interpolation
    Input data must be a list of lists of LDResult objects
    Example: [L1, L2, ....]
             where L1 = [LDResult1, LDResult2, ...]
                   L2 = [LDResult1, LDResult2, ...]
             Each item of L1 (or L2 ...) is a LDResult object
'''
def precision_recall_threshold(data, ldcut = 1):
    # Number of threshold points = number of SNPs in each list
    nsnps = max([len(x) for x in data])
    xvals = np.linspace(0, nsnps, nsnps)
    y1 = list()
    y2 = list()

    for d in data:
        _, _tpr, _ppv, _nsel, _, _, _, _ = roc.confusion_matrix(d, ldcut)
        tpr_est = np.interp(xvals, _nsel, _tpr)
        ppv_est = np.interp(xvals, _nsel, _ppv)
        y1.append(tpr_est)
        y2.append(ppv_est)

    recall        = np.array(y1).mean(axis = 0)
    precision     = np.array(y2).mean(axis = 0)
    recall_err    = np.array(y1).std(axis = 0)
    precision_err = np.array(y2).std(axis = 0)

    return xvals, precision, recall, precision_err, recall_err
