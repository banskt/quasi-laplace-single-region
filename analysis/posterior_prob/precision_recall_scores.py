import numpy as np

def confusion_matrix(ldresult, ldcut = 1.0):
    nitems = len(ldresult)
    ypred = np.array([x.stat for x in ldresult])
    ytrue = np.array([x.causality for x in ldresult])
    isort = np.argsort(ypred)[::-1]
    
    pos = np.sum(ytrue)
    neg = nitems - pos
    tp = 0	# True positives
    fp = 0	# False positives
    tpld = 0	# True positives (if in strong LD)
    fpld = 0    # False positives (not in LD with actual TP)

    tplist = list()
    fplist = list()
    tpldlist = list()
    fpldlist = list()
    alpha = np.max(ypred) + 1.0 # set threshold above the maximum value of predictions
    
    for j in range(nitems):
        if not ypred[isort[j]] == alpha:
            tplist.append(tp)
            fplist.append(fp)
            tpldlist.append(tpld)
            fpldlist.append(fpld)
            alpha = ypred[isort[j]]
        if ytrue[isort[j]] == 1:
            tp += 1
        else:
            fp += 1
        if ldresult[isort[j]].ld > ldcut:
            tpld += 1
        else:
            fpld += 1

    fplist.append(fp)
    tplist.append(tp)
    tpldlist.append(tpld)
    fpldlist.append(fpld)

    tpr = np.array([x / pos for x in tplist]) 							# TPR = Recall = TP / Negatives
    fpr = np.array([x / neg for x in fplist]) 							# FPR = FP / Negatives
    ppv = np.array([x[0] / sum(x) if sum(x) > 0 else 1 for x in zip(tplist, fplist)]) 		# PPV = Precision = TP / (TP + FP)
    fdr = np.array([x[1] / sum(x) if sum(x) > 0 else 0 for x in zip(tplist, fplist)]) 		# FDR = FP / (TP + FP)
    nsel = np.array([sum(x) for x in zip(tplist, fplist)]) 					# Number of y selected at each threshold
    ldtpr = np.array([x / pos for x in tpldlist])
    ldppv = np.array([x[0] / sum(x) if sum(x) > 0 else 0 for x in zip(tpldlist, fpldlist)])
    ldnsel = np.array([sum(x) for x in zip(tpldlist, fpldlist)])

    return fpr, tpr, ppv, nsel, fdr, ldtpr, ldppv, ldnsel
