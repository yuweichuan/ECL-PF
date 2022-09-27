import csv
import numpy as np
import pickle
import os
import time

# def protein_score_filtering(result, proportion=0.1):
#     """input is multiprocessing result from concurrent feature,
#     we suppose the larger the score, the more confident the spectrum,
#     output is the unique ID for each spectrum"""
#     proteinMap = {}
#     extraction = []  # extract the needed element
#     result = list(result)
#     for i in result:  # construct the protein score map
#         for j in i:
#             _topAlpha = j[0]
#             _topBeta = j[1]
#             for alpha in _topAlpha:
#                 extraction.append((alpha[0], alpha[3]))  # (score, description)
#             for beta in _topBeta:
#                 extraction.append((beta[0], beta[3]))  # (score, description)
#     extraction.sort(reverse=True)
#     length = int(len(extraction) * 0.1)
#     for i in range(length):  # use top scores to construct the protein map
#         proteins = extraction[i][1].split('$')
#         for protein in proteins:
#             proteinMap[protein] = proteinMap.setdefault(protein, 0) + extraction[i][0]
#     '''pick the ID from the top list'''
#
#     res = []
#
#     for i in result:
#         for j in i:
#             alphaList = []
#             betaList = []
#             a2 = []  # to calculate average score of alpha candidates
#             b2 = []  # to calculate average score of beta candidates
#             aCurrMax = -1E12  # initial  max score
#             bCurrMax = -1E12  # initial  max score
#             for alpha in j[0]:
#                 a2.append(alpha[0])
#                 proteins = alpha[3].split('$')
#                 pScore = sum([proteinMap.setdefault(i, 0) for i in proteins]) / len(proteins)
#                 if pScore > aCurrMax:
#                     aCurrMax = pScore
#                     alphaList = [alpha]
#                 elif pScore == aCurrMax:
#                     alphaList.append(alpha)
#             for beta in j[1]:
#                 b2.append(beta[0])
#                 proteins = beta[3].split('$')
#                 pScore = sum([proteinMap.setdefault(i, 0) for i in proteins]) / len(proteins)
#                 if pScore > bCurrMax:
#                     bCurrMax = pScore
#                     betaList = [beta]
#                 elif pScore == bCurrMax:
#                     betaList.append(beta)
#             a2 = sum(a2) / len(a2)
#             b2 = sum(b2) / len(b2)
#             for alpha in alphaList:
#                 for beta in betaList:
#                     res.append({'CID_scan': j[3], 'mass': j[4], 'marker': j[2], 'alpha1': list(alpha),
#                                 'beta1': list(beta), 'alpha2': a2, 'beta2': b2})
#     return res


def truncate(ls1, k):
    res = sorted(ls1)
    s1 = max(-k, -len(ls1))
    while s1 > -len(ls1):
        if ls1[s1][0] == ls1[s1 - 1][0]:
            s1 -= 1
        else:
            break
    return res[s1:]


def elbow(_extr, pt=0.01):
    """determine the ReLu shape elbow, require sorted from largest to smallest"""
    # _extr = sorted(extr, key=lambda x: -x[1][2])  # largest to smallest
    T = 0
    D = 0
    pepFdr = []
    for _i in _extr:
        if any(np.array(_i[1][1].split('$')).astype('int') % 2):
            T += 1
        else:
            D += 1
        if T == 0:
            pepFdr.append(1)
        else:
            pepFdr.append(D / T)
    tt = [ii[1][2] for ii in _extr if any(np.array(ii[1][1].split('$')).astype('int') % 2)]
    dd = [ii[1][2] for ii in _extr if not any(np.array(ii[1][1].split('$')).astype('int') % 2)]

    for _i in range(len(pepFdr)-1, -1, -1):
        if pepFdr[_i] <= pt:
            plt.hist(tt, alpha=0.5)
            plt.hist(dd, alpha=0.5)
            plt.scatter(_extr[_i][1][2], 0)
            plt.show()
            return _extr[_i][1][2]
    return _extr[0][1][2]


def skewed_tanh(x, k):
    x = np.array(x)
    x[x > 400] = 400  # when k = 1, x = 400 almost reaches the upperbound 99.93, this is to avoid overflow
    a = 100 / k
    b = np.log((100 - k) / (100 + k)) + 1
    return a * (np.exp(x) - np.exp(b * x)) / (np.exp(x) + np.exp(b * x))


def merge_seq(ls):
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
    """input is multiprocessing result from concurrent feature,
    we suppose the larger the score, the more confident the spectrum,
    output is the unique ID for each spectrum"""
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
    # elb = 2 * np.median([s[1][2] for s in extraction])  # two-times median number
    elb_ls = [s[1][2] for s in extraction]
    if not elb_ls:
        elb = 0
    else:
        elb = np.percentile(elb_ls, 98)  # top 2% number

    # elb = elbow(extraction)
    # length = int(len(extraction) * proportion)

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


def export_all(result):
    """input is multiprocessing result from concurrent feature,
    export everything first, then to determine the result"""
    proteinMap = {}
    extraction = []  # extract the needed element
    result = list(result)
    res = []
    for i in result:  # construct the protein score map
        for j in i:
            cur = []
            a2 = []  # to calculate average score of alpha candidates
            b2 = []  # to calculate average score of beta candidates
            _topAlpha = j[0]
            _topBeta = j[1]
            for alpha in _topAlpha:
                a2.append(alpha[0])
                for beta in _topBeta:
                    b2.append(beta[0])
                    cur.append({'CID_scan': j[3], 'mass': j[4], 'marker': j[2], 'alpha1': list(alpha),
                                'beta1': list(beta)})
            a2 = sum(a2) / len(a2)
            b2 = sum(b2) / len(b2)
            for ele_cur in cur:
                ele_cur['alpha2'] = a2
                ele_cur['beta2'] = b2
            res.extend(cur)
    return res


def raw_file2_unique(path, proportion=0.1):
    proteinMap = {}
    extraction = set()
    scanDict = {}
    with open(path) as csvfile:
        file = csv.DictReader(csvfile)
        proA = set()
        proB = set()
        for row in file:

            if scanDict.setdefault(row['CID_scan'], []) == list():
                scanDict[row['CID_scan']] = [{(row['alpha'], row['a_protein'], row['a1_score'])},
                                             {(row['beta'], row['b_protein'], row['b1_score'])}]
            else:
                scanDict[row['CID_scan']][0].add((row['alpha'], row['a_protein'], row['a1_score']))
                scanDict[row['CID_scan']][1].add((row['beta'], row['b_protein'], row['b1_score']))
            desA = row['a_protein'].split(';')
            desB = row['b_protein'].split(';')
            for des in desA:
                if (des, row['CID_scan']) not in proA:
                    proA.add((des, row['CID_scan']))
                    extraction.add((float(row['a1_score']), des, row['CID_scan']))
            for des in desB:
                if (des, row['CID_scan']) not in proB:
                    proB.add((des, row['CID_scan']))
                    extraction.add((float(row['b1_score']), des, row['CID_scan']))
        extraction = list(extraction)
        extraction.sort(reverse=True)
        length = int(len(extraction) * proportion)
        for i in range(length):  # use top scores to construct the protein map
            protein = extraction[i][1]
            proteinMap[protein] = proteinMap.setdefault(protein, 0) + extraction[i][0]
    # plt.hist(proteinMap.values(), bins=100)
    # plt.show()
    res = []
    proteinMap['asb'] = -1
    for scan, ab in scanDict.items():
        a1 = ('pep','asb')
        a1s = -1
        b1 = ('pep', 'asb')
        b1s = -1
        for alpha in ab[0]:
            ps = alpha[1].split(';')
            ps = [proteinMap.setdefault(ii, 0) for ii in ps]
            ps = sum(ps) / len(ps)
            if ps > a1s:
                a1s = ps
                a1 = [alpha]
            elif ps == a1s:
                a1.append(alpha)
        for beta in ab[1]:
            ps = beta[1].split(';')
            ps = [proteinMap.setdefault(ii, 0) for ii in ps]
            ps = sum(ps) / len(ps)
            if ps > b1s:
                b1s = ps
                b1 = [beta]
            elif ps == b1s:
                b1.append(beta)
        # if len(a1)*len(b1) < 60:
        for aa in a1:
            for bb in b1:
                res.append((scan, aa, bb))
    with open(r"D:\OneDrive - HKUST Connect\Research\Reading Report\20210910\simulated dataset_t10_t10_noEtdScore\extractAll\Xcorr\test.csv", 'w', newline='') as csvfile:
        fieldnames = ['CID_scan', 'alpha', 'beta', 'a_protein', 'b_protein', 'a1_score', 'b1_score', 're_rank_score',
                      'mass', 'pos_a', 'pos_b', 'marker', 'a2_score', 'b2_score', 'a_mass', 'b_mass']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for spectrum in res:

            writer.writerow(
                {'CID_scan': spectrum[0], 'alpha': spectrum[1][0], 'a_protein': spectrum[1][1], 'a1_score': spectrum[1][2],
                 'beta': spectrum[2][0], 'b_protein': spectrum[2][1], 'b1_score': spectrum[2][2],
                 're_rank_score': float(spectrum[1][2]) * float(spectrum[2][2]), 'mass':1000, 'pos_a':1, 'pos_b':1,
                 'marker':4, 'a2_score':1, 'b2_score':1, 'a_mass':1, 'b_mass':1})


if __name__ == '__main__':
    '''determine how the protein score functions'''
    '''extract results'''
    # ls = dict()
    # path = r"C:\Users\czhouau\Desktop\MS181452_F18_HCD.csv"
    # with open(path, newline='') as csvfile:
    #     reader = csv.DictReader(csvfile)
    #     for row in reader:
    #         ls.setdefault(int(row['CID_scan']), [set(), set()])[0].add((round(float(row['a1_score']), 3), row['alpha'], row['a_protein']))
    #         ls.setdefault(int(row['CID_scan']), [set(), set()])[1].add((round(float(row['b1_score']), 3), row['beta'], row['b_protein']))
    #
    # print('done export')
    # with open('proteinscoreDB', 'wb') as f:
    #     pickle.dump(ls, f)
    '''label ID'''
    # proteins = set()
    # with open(r"C:\Users\czhouau\Desktop\MS181452_F18_HCD_inter.csv") as f:
    #     a = csv.DictReader(f)
    #     for row in a:
    #         alpha = row['a_protein'].split('; ')
    #         beta = row['b_protein'].split('; ')
    #         for p in alpha+beta:
    #             proteins.add(p)
    #
    # with open(r"C:\Users\czhouau\Desktop\MS181452_F18_HCD_intra.csv") as f:
    #     a = csv.DictReader(f)
    #     for row in a:
    #         alpha = row['a_protein'].split('; ')
    #         beta = row['b_protein'].split('; ')
    #         for p in alpha+beta:
    #             proteins.add(p)
    # with open('labelID','wb') as f:
    #     pickle.dump(proteins, f)
    '''solid prove of PFM'''
    # with open('labelID', 'rb') as f:
    #     pid = pickle.load(f)
    # with open('proteinscoreDB', 'rb') as f:
    #     ls = pickle.load(f)
    # print(len(ls))
    # print(ls[10898][1])
    # print(len(ls[10898][1]))
    # tp = set()
    # T = dict()
    # D = dict()
    # for ele in ls.values():
    #     ele0 = max(ele[0])
    #     ele1 = max(ele[1])
    #     # if ele0[1] not in tp:
    #     #     tp.add(ele0[1])
    #     for i in ele0[2].split(';'):
    #         if i not in pid:
    #             D[i] = D.setdefault(i, 0) + 1
    #         else:
    #             T[i] = T.setdefault(i, 0) + 1
    #     # if ele1[1] not in tp:
    #     #     tp.add(ele1[1])
    #     for i in ele1[2].split(';'):
    #         if i not in pid:
    #             D[i] = D.setdefault(i, 0) + 1
    #         else:
    #             T[i] = T.setdefault(i, 0) + 1
    # print(np.mean(list(T.values())))
    # print(np.mean(list(D.values())))
    # sns.distplot(list(T.values()))
    # sns.distplot(list(D.values()))
    # # plt.hist(D, [1 for i in range(len(D))])
    # plt.show()
    # print(len(T), len(D))

    '''PFM calculation'''
    pathList = [r"C:\Users\czhouau\Desktop\New folder\F1_Agilent5_20160430_DF_InsolTh22_alignment"]
    ms1 = 1e-5
    xlMass = 158.0038
    for path in pathList:
        with open(path, 'rb') as f:
            ls = pickle.load(f)  # (alphas, betas, #marker, scan, preMass)
        with open('protein_name_heck', 'rb') as f:
            pname = pickle.load(f)

        maxls = []
        path = path.split('\\')[-1]
        for ele in ls:

            max1 = max(ele[0])[0]
            max2 = max(ele[1])[0]
            ele0 = [i for i in ele[0] if i[0] == max1]
            ele1 = [i for i in ele[1] if i[0] == max2]
            re0 = len(ele0)
            re1 = len(ele1)
            ele0 = [(re0, i) for i in ele0]
            ele1 = [(re1, i) for i in ele1]
            # if max1 > 0:
            if max1 > 0:
                maxls.extend(ele0)  # (#tying, (score, seq, pos, #des$#des, chainMass))
            # if max2 > 0:
            if max2 > 0:
                maxls.extend(ele1)
        '''protein database: protein with identified peptides'''
        # elb = 2 * np.median([ii[1][0] for ii in maxls])
        elb = np.percentile([ii[1][0] for ii in maxls], 98)  # top 10% number
        maxls = [ii for ii in maxls if ii[1][0] >= elb]
        pdb = dict()
        for ele in maxls:
            proteins = list(map(int, ele[1][3].split('$')))  # in integer format
            for protein in proteins:  # need norm/#tying
                pdb[protein] = pdb.setdefault(protein, [])
                pdb[protein].append((1/ele[0], ele[1][1]))

        for i in pdb.keys():
            pdb[i] = sum(skewed_tanh(merge_seq(pdb[i]), pname[i][1]))
            # pdb[i] = [sum(skewed_tanh(merge_seq(pdb[i]))), pdb[i]]
        # T = []
        # D = []
        # for i,j in pdb.items():
        #     if 'DECOY' in i:
        #         D.append([i,j])
        #     else:
        #         T.append([i,j])
        # D = sorted(D, key=lambda x: (-x[1][0]))
        # T = sorted(T, key=lambda x: (-x[1][0]))
        # for i in range(100):
        #     print(D[i])
        # with open('tanh', 'wb') as f:
        #     pickle.dump(pdb, f)
        pts = []
        tops = [1]
        proportion = [1]
        # tops = [30]
        # proportion = [0.4]


        # tops = [1, 5, 10, 15]
        # proportion = [0, 0.5, 0.9, 1]
        for top in tops:
            for prop in proportion:
                result = []
                for i in ls:
                    cura = -10
                    curb = -10

                    max0 = max(i[0])[0]
                    max1 = max(i[1])[0]
                    # if max0 <= 0 and max1 <= 0:
                    if max0 < 0 or max1 < 0:  # -20, -20 for merox; 0, 0 for XlinkX and Xcorr
                        pass
                    else:
                        # for alpha in truncate([ii for ii in i[0] if ii[0] >= prop * max0], top):
                        for alpha in truncate([ii for ii in i[0] if ii[0] >= prop * max0], top):
                            proteins = list(map(int, alpha[3].split('$')))
                            # pScore = sum([pdb.setdefault(ii, 0) for ii in proteins]) / len(proteins)
                            pScore = max([pdb.setdefault(ii, 0) for ii in proteins])
                            if pScore > cura:
                                cura = pScore
                                alphaList = [alpha]
                            elif pScore == cura:
                                alphaList.append(alpha)

                        # for beta in truncate([ii for ii in i[1] if ii[0] >= prop * max1], top):
                        for beta in truncate([ii for ii in i[1] if ii[0] >= prop * max1], top):
                            proteins = list(map(int, beta[3].split('$')))
                            # pScore = sum([pdb.setdefault(ii, 0) for ii in proteins]) / len(proteins)
                            pScore = max([pdb.setdefault(ii, 0) for ii in proteins])
                            if pScore > curb:
                                curb = pScore
                                betaList = [beta]
                            elif pScore == curb:
                                betaList.append(beta)

                        result.append((i[2], i[3], i[4], alphaList, betaList, cura, curb))  # marker, scan, preMass, alphas, betas, pa, pb

                with open('{}.csv'.format(path), 'w', newline='') as csvfile:
                    fieldnames = ['CID_scan', 'mass', 'marker', 're_rank_score', 'alpha', 'pos_a', 'a1_score', 'a2_score',
                                  'a_mass', 'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein',
                                  'q_value']
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

                    writer.writeheader()
                    for spectrum in result:
                        for alpha in spectrum[3]:
                            for beta in spectrum[4]:
                                if np.isclose(spectrum[2], alpha[4] + beta[4] + xlMass, rtol=ms1):
                                    a_p = [pname[ii][0].split(' ')[0] for ii in list(map(int, alpha[3].split('$')))]
                                    a_p = [ii for ii in a_p if 'DECOY' not in ii] if [ii for ii in a_p if
                                                                                      'DECOY' not in ii] else a_p
                                    a_p = ';'.join(a_p)
                                    b_p = [pname[ii][0].split(' ')[0] for ii in list(map(int, beta[3].split('$')))]
                                    b_p = [ii for ii in b_p if 'DECOY' not in ii] if [ii for ii in b_p if
                                                                                      'DECOY' not in ii] else b_p
                                    b_p = ';'.join(b_p)
                                    writer.writerow(
                                        {'CID_scan': int(spectrum[1]), 'mass': spectrum[2], 'marker': spectrum[0],
                                         're_rank_score': (spectrum[5] * spectrum[6]), 'alpha': alpha[1],
                                         'pos_a': alpha[2],
                                         'a1_score': alpha[0], 'a_mass': alpha[4],
                                         'a_protein': a_p,
                                         'beta': beta[1], 'pos_b': beta[2], 'b1_score': beta[0],  #beta[0]
                                         'b_protein': b_p,
                                         'q_value': 1, 'b_mass': beta[4],
                                         'a2_score': spectrum[5], 'b2_score': spectrum[6]})
                os.system('python splitctrl.py {}.csv'.format(path))
                time.sleep(2)
                with open('{}_final.csv'.format(path)) as f:
                    intra = csv.DictReader(f)
                    intraLen = len(list(intra))


                print('top is {}, proportion is {}, number is {}'.format(top, prop, intraLen))
                pts.append([top, prop, intraLen])

        pts = np.array(pts)
        with open(path,'wb') as f:
            pickle.dump(pts, f)
        index = np.argmax(pts[:, 2])
        print(pts[index])

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2])
        ax.set_xlabel('#top')
        ax.set_ylabel('Proportion')
        ax.set_zlabel('#CSM')
        plt.savefig('{}.png'.format(path))
        plt.close(fig)
