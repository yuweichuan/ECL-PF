import sys
import csv
import numpy as np
from pyteomics import fasta
import os


def separate(data):
    """
    Separate the ran out data as inter and intra-protein groups.
    Inter group includes homo-oligomers
    :param data:
    """
    inter = []  # initialize inter group
    intra = []  # initialize intra group
    for sub_data in data:
        if set(sub_data['a_protein']) & set(sub_data['b_protein']):
            sub_data['a_protein'] = list(set(sub_data['a_protein']) & set(sub_data['b_protein']))
            sub_data['b_protein'] = list(set(sub_data['a_protein']) & set(sub_data['b_protein']))
            intra.append(dict(sub_data))
        else:
            if len(sub_data['a_protein']) == 1 and sub_data['a_protein'] == sub_data['b_protein']:
                pass
            else:
                inter.append(dict(sub_data))
    return inter, intra


def inter_fdr(inter_data, fdr=0.01):
    """
    Calculate the FDR by the separate FDR control method,
    the function is adapted to both inter and intra-protein groups
    :param inter_data:
    :param fdr:
    """
    data = sorted(inter_data, key=lambda x: -(x['b1_score'] + x['a1_score']))  # sort the CSMs by the descending order
    TT = 0  # target-target set
    TD = 0  # target-decoy set
    DD = 0  # decoy-decoy set
    for sub in data:
        tga = True if 'DECOY' in sub['a_protein'][0] else False
        tgb = True if 'DECOY' in sub['b_protein'][0] else False
        if tgb and tga:
            DD += 1 / weightList[sub['CID_scan']]
        elif (not tga) and (not tgb):
            TT += 1 / weightList[sub['CID_scan']]
        else:
            TD += 1 / weightList[sub['CID_scan']]
        if TT == 0:
            sub['fdr'] = 1
        else:
            if TD > DD:
                sub['fdr'] = (TD - DD) / TT
            else:
                sub['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(data) - 1, -1, -1):
        if data[idx]['fdr'] < minfdr:
            minfdr = data[idx]['fdr']
        data[idx]['q_value'] = minfdr

    data = [sub_data for sub_data in data if sub_data['q_value'] <= fdr and 'DECOY' not in sub_data['a_protein'][0] and
            'DECOY' not in sub_data['b_protein'][0]]
    return data


def tanhSkew(x, k):
    """
    Counting function
    :param x:
    :param k:
    """
    base = (100.0 - k) / (100.0 + k)
    p = base ** x
    res = (100.0 / k) * (1.0 - p) / (1.0 + p)
    return res


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


class Record:
    def __init__(self, scan, mass, score, peptide1, xlSite1, pepMass1, protein1, score1, peptide2, xlSite2, pepMass2,
                 protein2, score2):
        self.scan = scan
        self.mass = mass
        self.score = score
        self.peptide1 = peptide1
        self.xlSite1 = xlSite1
        self.pepMass1 = pepMass1
        self.protein1 = protein1
        self.score1 = score1
        self.peptide2 = peptide2
        self.xlSite2 = xlSite2
        self.pepMass2 = pepMass2
        self.protein2 = protein2
        self.score2 = score2


if __name__ == '__main__':
    records = dict()
    csmPath = sys.argv[1]
    check_file = os.path.exists(csmPath)
    if not check_file:
        quit()
    proteinPath = sys.argv[2]
    FDR = float(sys.argv[3])

    '''Necessary steps in the protein feedback, construct numberK'''
    protein_name = dict()
    upperScore = 100  # upper limit protein score
    for des, seq in fasta.read(proteinPath):
        key = des.split(' ')[0]
        dkey = 'DECOY_' + key
        value = min(upperScore - 1, seq.count('K') + 1)
        protein_name[key] = value
        protein_name[dkey] = value

    '''Construct the output results'''
    csv.field_size_limit(500 * 1024 * 1024)
    with open(csmPath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            key = int(row['ScanNum']) * 1E4 + float(row['Mass'])
            cur = Record(int(row['ScanNum']), float(row['Mass']), float(row['Score']),
                         row['Peptide#1'].split('.')[1], int(row['LinkSite#1']), float(row['PepMass#1']), row['Protein#1'].split('; '), float(row['Score#1']),
                         row['Peptide#2'].split('.')[1], int(row['LinkSite#2']), float(row['PepMass#2']), row['Protein#2'].split('; '), float(row['Score#2']))
            records.setdefault(key, []).append(cur)

    '''Construct protein score database'''
    extraction = []  # used to determine the high confident results in the histogram
    for value in records.values():
        curMax = value[-1].score
        temp = []  # [(peptide, [des1 des2 des3], score), (), (), ...]
        tieTo1 = 0
        for ele in value[::-1]:
            if ele.score < curMax:
                break
            else:
                tieTo1 += 1
                seq1 = ''.join([i for i in ele.peptide1 if i.isupper()])
                seq2 = ''.join([i for i in ele.peptide2 if i.isupper()])
                temp.append((seq1, ele.protein1, ele.score1))
                temp.append((seq2, ele.protein2, ele.score2))

        for ele in temp:
            extraction.append((ele[2], tieTo1, ele[1], ele[0]))  # [(pep score, tieTo1, des des des, pep seq), (),...]

    '''Find score cutoff for the high confident peptides'''
    elb_ls = [s[0] for s in extraction]
    maxNumber = 0
    final_output = []
    for percent in np.arange(98, 98.4, 0.5):  # could be used to iteratively determine the best results, currently with 98 only
        if not elb_ls:
            elb = 0
        else:
            elb = np.percentile(elb_ls, percent)  # top 2% number

        '''Construct protein score database'''
        proteinMap = dict()
        proteinScore = dict()
        for ele in extraction:
            for pName in ele[2]:
                proteinMap.setdefault(pName, []).append(((1 if ele[0] >= elb else 0) / ele[1], ele[3]))  # (norm/#tying, seq)

        for key in proteinMap.keys():
            pl = merge_seq(proteinMap[key])
            ps = 0
            for p in pl:
                ps += tanhSkew(p, protein_name[key])
            proteinScore[key] = ps

        '''pick the right count'''
        dataAll = []
        for ele in records.values():
            maxPS = 1
            curVec = []
            for record in ele[::-1]:
                alphaPS = max([proteinScore.setdefault(name, 0) for name in record.protein1])
                betaPS = max([proteinScore.setdefault(name, 0) for name in record.protein2])
                if maxPS < (alphaPS + 1.0) * (betaPS + 1.0):
                    curVec = [record]
                    maxPS = (alphaPS + 1.0) * (betaPS + 1.0)
                elif maxPS == (alphaPS + 1.0) * (betaPS + 1.0):
                    curVec.append(record)
            '''only retrieve the highest scored ones'''
            maxScore = curVec[0].score
            for cur in curVec:
                if maxScore <= cur.score:
                    element = {'CID_scan': cur.scan, 'Mass': cur.mass, 'Score': cur.score, 'PF_Score': maxPS,
                               'alpha': cur.peptide1, 'pos_a': cur.xlSite1, 'a1_score': cur.score1, 'a_protein': cur.protein1,
                               'a_mass': cur.pepMass1, 'b_mass': cur.pepMass2, 'beta': cur.peptide2, 'pos_b': cur.xlSite2,
                               'b1_score': cur.score2, 'b_protein': cur.protein2}
                    dataAll.append(element)
                else:
                    break

        'resort by scan number'
        dataAll = sorted(dataAll, key=lambda x: (x['CID_scan']))
        '''remove redundant decoys'''
        cur_scan = -1
        currD = []
        D = []
        redundantFlag = False
        for idx in range(len(dataAll)):
            if dataAll[idx]['CID_scan'] != cur_scan:
                '''add last to D'''
                if redundantFlag:
                    D.extend(currD)
                '''initiate'''
                cur_scan = dataAll[idx]['CID_scan']
                currD = []
                redundantFlag = False
                if 'DECOY' in dataAll[idx]['a_protein'][0] or 'DECOY' in dataAll[idx]['b_protein'][0]:
                    currD.append(idx)
            else:
                redundantFlag = True
                if 'DECOY' in dataAll[idx]['a_protein'][0] or 'DECOY' in dataAll[idx]['b_protein'][0]:
                    currD.append(idx)
        dataAll = [dataAll[ii] for ii in range(len(dataAll)) if ii not in D]

        '''distinguish ambiguity by related ones if we can'''
        solid = []  # unique results
        redundant = []  # redundant results
        lastScan = -1
        popFlag = False
        for ele in dataAll:
            if ele['CID_scan'] != lastScan:
                lastScan = ele['CID_scan']
                popFlag = False
                solid.append(ele)
            else:
                redundant.append(ele)
                if not popFlag:
                    redundant.append(solid.pop())
                    popFlag = True
        solidSet = set()

        for ele in solid:
            pep1 = ''.join([ii for ii in ele['alpha']if ii.isupper()])
            pep2 = ''.join([ii for ii in ele['beta']if ii.isupper()])
            solidSet.add((pep1, pep2))
            solidSet.add((pep2, pep1))
        removeFlag = False
        cur_scan = -1
        removeList = []
        curList = []
        for ele in redundant:
            if cur_scan != ele['CID_scan']:
                cur_scan = ele['CID_scan']
                if removeFlag:
                    removeList.extend(curList)
                curList = []

            pep1 = ''.join([ii for ii in ele['alpha'] if ii.isupper()])
            pep2 = ''.join([ii for ii in ele['beta'] if ii.isupper()])
            if (pep1, pep2) in solidSet:
                removeFlag = True
            else:
                curList.append(ele)
        if removeFlag:
            removeList.extend(curList)
        redundant = [i for i in redundant if i not in removeList]
        dataAll = redundant + solid

        '''this is to distribute the weights of redundant scans'''
        weightList = {}
        for ele in dataAll:
            weightList[ele['CID_scan']] = weightList.setdefault(ele['CID_scan'], 0) + 1

        '''divide each type into intra and inter'''
        inter, intra = separate(dataAll)
        res_inter = inter_fdr(inter, fdr=FDR)
        res_intra = inter_fdr(intra, fdr=FDR)

        '''get rid of sequence identified only once with 0 score in intra'''

        seqDict = dict()
        eliminate = []
        for ele in res_intra:
            seqDict.setdefault(ele['alpha'], list()).append((ele['CID_scan'], ele['a1_score']))
            seqDict.setdefault(ele['beta'], list()).append((ele['CID_scan'], ele['b1_score']))
        for ele in seqDict.values():
            if len(ele) == 1 and ele[0][1] == 0:
                eliminate.append(ele[0][0])

        '''store by one final csv file'''
        final_res = []
        inter_scan = dict()
        res_index = 0

        for ele in res_inter:
            ele['type'] = 'inter'
            final_res.append(ele)
            inter_scan[ele['CID_scan']] = res_index
            res_index += 1
        for ele in res_intra:
            if ele['CID_scan'] in eliminate:
                continue
            if ele['CID_scan'] in inter_scan.keys():
                if len(final_res[inter_scan[ele['CID_scan']]]['a_protein']) == 1 and final_res[inter_scan[ele['CID_scan']]]['a_protein'] == final_res[inter_scan[ele['CID_scan']]]['b_protein']:

                    final_res[inter_scan[ele['CID_scan']]]['type'] = r'intra'
                    final_res[inter_scan[ele['CID_scan']]]['q_value'] = min(ele['q_value'],
                                                                            final_res[inter_scan[ele['CID_scan']]]['q_value'])
                else:
                    final_res[inter_scan[ele['CID_scan']]]['type'] = r'intra/inter'
                    final_res[inter_scan[ele['CID_scan']]]['q_value'] = min(ele['q_value'],
                                                                            final_res[inter_scan[ele['CID_scan']]]['q_value'])
            else:
                ele['type'] = 'intra'
                final_res.append(ele)

            '''calculate average number of one cl peptides being identified, remove less confident by '''
        new_res = []
        indivNum = dict()
        pepNum = dict()
        for ele in final_res:
            pepNum['{}-{}'.format(ele['alpha'], ele['beta'])] = \
                pepNum.setdefault('{}-{}'.format(ele['alpha'], ele['beta']), 0) + 1

        if list(pepNum.values()):
            aveNum = np.mean(list(pepNum.values()))
        else:
            aveNum = 0

        for ele in final_res:

            if pepNum['{}-{}'.format(ele['alpha'], ele['beta'])] >= aveNum / 4:
                new_res.append(ele)
        final_res = sorted(new_res, key=lambda x: -x['b1_score'] - x['a1_score'])
        
        if maxNumber < len(final_res):
            maxNumber = len(final_res)
            final_output = final_res

    '''write the output file'''
    with open("{}_final.csv".format(csmPath.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
        fieldnames = ['CID_scan', 'Mass', 'Score', 'PF_Score', 'alpha', 'pos_a', 'a_mass', 'a1_score', 'a_protein', 'beta', 'pos_b', 'b_mass', 'b1_score', 'b_protein', 'type', 'q_value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for spectrum in final_output:
            a_protein = '; '.join([i.split(' ')[0] for i in set(spectrum['a_protein'])])
            b_protein = '; '.join([i.split(' ')[0] for i in set(spectrum['b_protein'])])
            writer.writerow(
                {'CID_scan': spectrum['CID_scan'], 'Mass': spectrum['Mass'], 'PF_Score': spectrum['PF_Score'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
                 'Score': spectrum['Score'], 'alpha': spectrum['alpha'],
                 'pos_a': spectrum['pos_a'], 'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
                 'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
                 'b_protein': b_protein, 'q_value': spectrum['q_value'],
                 'type': spectrum['type']})
    os.remove(csmPath)
