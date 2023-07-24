import numpy as np


def skewed_tanh(x, k):
    """
    Counting function
    :param x:
    :param k:
    """
    x = np.array(x)
    x[x > 400] = 400  # when k = 1, x = 400 almost reaches the upperbound 99.93, this is to avoid overflow
    a = 100 / k
    b = np.log((100 - k) / (100 + k)) + 1
    return a * (np.exp(x) - np.exp(b * x)) / (np.exp(x) + np.exp(b * x))


def merge_seq(ls):
    """
    Calculate the protein score
    :param ls:
    """
    res = []
    seqls = []
    for _l in ls:
        if _l[1] not in seqls:
            seqls.append(_l[1])
            res.append(_l[0])
        else:
            idx = seqls.index(_l[1])
            res[idx] += _l[0]
    return res


def protein_score_filtering(_result, tol1, xlMass, numberK):
    """
    Protein feedback process. input is multiprocessing result from concurrent feature,
    We suppose the larger the score, the more confident the spectrum,
    output is the unique ID for each spectrum
    :param _result:
    :param tol1:
    :param xlMass:
    :param numberK:
    """
    proteinMap = {}
    extraction = []  # extract the needed element
    _result = list(_result)
    pickle_res = []

    for i in _result:  # construct the protein score map
        for j in i:
            pickle_res.append(j)
            _topAlpha = j[0]
            _topBeta = j[1]

            maxAlpha = max(_topAlpha)[0]
            if maxAlpha > 0:
                _alpha = [(ii[1], ii[3], ii[0]) for ii in _topAlpha if ii[0] == maxAlpha]  # (seq, des$des$des, score)
                lenA = len(_alpha)
                _alpha = [(lenA, ii) for ii in _alpha]  # (#tying the 1st, (seq, des$des$des, score))
                extraction.extend(_alpha)

            maxBeta = max(_topBeta)[0]
            if maxBeta > 0:
                _beta = [(ii[1], ii[3], ii[0]) for ii in _topBeta if ii[0] == maxBeta]  # (seq, des$des$des, score)
                lenB = len(_beta)
                _beta = [(lenB, ii) for ii in _beta]
                extraction.extend(_beta)

    extraction = sorted(extraction, key=lambda x: -x[1][2])
    if len(extraction) > 0:
        elb = np.percentile([s[1][2] for s in extraction], 98)  # top 2% number
    else:
        elb = 100

    for i in range(len(extraction)):  # use top scores to construct the protein map

        proteins = extraction[i][1][1].split('$')
        for protein in proteins:
            proteinMap.setdefault(int(protein), []).append(((1 if extraction[i][1][2] / elb >= 1 else 0) / extraction[i][0],
                                                            extraction[i][1][0]))  # (norm/#tying, seq)
            # proteinMap.setdefault(int(protein), []).append((0, extraction[i][1][0]))  # (norm/#tying, seq)
    '''pick the ID from the top list'''
    for i in proteinMap.keys():
        proteinMap[i] = sum(skewed_tanh(merge_seq(proteinMap[i]), numberK[i]))
    res = []
    for i in _result:
        for j in i:
            alphaList = []
            betaList = []
            a2 = []  # to calculate average score of alpha candidates
            b2 = []  # to calculate average score of beta candidates
            aCurrMax = -10  # initial  max score -10, synthetic 1e-3
            bCurrMax = -10  # initial  max score -10, synthetic 1e-3

            '''the following four lines, maxA and maxB are important to reduce computing time'''
            maxA = max(j[0])[0]
            maxB = max(j[1])[0]
            if maxA < 0 or maxB < 0:
                continue

            for alpha in j[0]:
                a2.append(alpha[0])
                proteins = alpha[3].split('$')
                pScore = max([proteinMap.setdefault(int(i), 0) for i in proteins])
                if pScore > aCurrMax:
                    aCurrMax = pScore
                    alphaList = [alpha]
                elif pScore == aCurrMax:
                    alphaList.append(alpha)
            for beta in j[1]:
                b2.append(beta[0])
                proteins = beta[3].split('$')
                pScore = max([proteinMap.setdefault(int(i), 0) for i in proteins])
                if pScore > bCurrMax:
                    bCurrMax = pScore
                    betaList = [beta]
                elif pScore == bCurrMax:
                    betaList.append(beta)
            a2 = sum(a2) / len(a2)
            b2 = sum(b2) / len(b2)
            for alpha in alphaList:
                for beta in betaList:
                    if np.isclose(j[4], xlMass + alpha[4] + beta[4], rtol=tol1, atol=0):
                        res.append({'protein_score': aCurrMax * bCurrMax, 'CID_scan': j[3], 'mass': j[4], 'marker': j[2], 'alpha1': list(alpha),
                                    'beta1': list(beta), 'alpha2': a2, 'beta2': b2})
    return res
