import sys
import csv
import numpy as np


def separate(data):
    """given the data with dict as ele, separate them into inter (stay original results) and intra"""
    inter = []
    intra = []
    for sub_data in data:
        inter.append(dict(sub_data))
        if set(sub_data['a_protein']) & set(sub_data['b_protein']):
            # sub_data['a_protein'] = list(set(sub_data['a_protein']) & set(sub_data['b_protein']))
            # sub_data['b_protein'] = list(set(sub_data['a_protein']) & set(sub_data['b_protein']))
            intra.append(dict(sub_data))
            
        # if len(sub_data['a_protein']) == 1 and sub_data['a_protein'] == sub_data['b_protein']:
        #     pass
        # else:
        #     inter.append(dict(sub_data))
    return inter, intra


def inter_fdr(inter_data, fdr=0.01):
    data = sorted(inter_data, key=lambda x: -(x['b1_score'] + x['a1_score']))
    
    TT = 0
    TD = 0
    DD = 0
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


if __name__ == '__main__':

    dataAll = []
    FDR = 0.01
    path = sys.argv[1]
    FDR = float(sys.argv[2])
    csv.field_size_limit(500 * 1024 * 1024)
    weightList = {}
    with open(path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ele = {'CID_scan': int(row['CID_scan']), 'mass': float(row['mass']),
                   'marker': float(row['marker']), 're_rank_score': float(row['re_rank_score']),
                   'alpha': row['alpha'], 'pos_a': int(row['pos_a']), 'a1_score': float(row['a1_score']),
                   'a2_score': float(row['a2_score']), 'a_mass': float(row['a_mass']),
                   'a_protein': row['a_protein'].split(';'),
                   'beta': row['beta'], 'pos_b': int(row['pos_b']), 'b1_score': float(row['b1_score']),
                   'b2_score': float(row['b2_score']), 'b_mass': float(row['b_mass']),
                   'b_protein': row['b_protein'].split(';')}
            dataAll.append(ele)

    'resort by scan number'
    dataAll = sorted(dataAll, key=lambda x: (-x['CID_scan'], -x['b1_score'] - x['a1_score']))

    '''first only retain the highest scored ones'''
    retain = []
    cur_scan = -1
    cur_max = -1E10
    curList = []
    for ele in dataAll:
        if ele['CID_scan'] != cur_scan:
            cur_max = ele['a1_score'] + ele['b1_score']
            retain.append(ele)
            cur_scan = ele['CID_scan']
            curList = [(ele['alpha'], ele['beta'])]
        else:
            if ele['a1_score'] + ele['b1_score'] == cur_max and (ele['alpha'], ele['beta']) not in curList:
                retain.append(ele)
                curList.append((ele['alpha'], ele['beta']))
    dataAll = retain

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
    solid = []
    redundant = []
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

    '''remove redundant (optional)'''
    redundantpool = {}
    for ele in dataAll:
        redundantpool[ele['CID_scan']] = redundantpool.setdefault(ele['CID_scan'], 0) + 1
    res_new = []

    for ele in dataAll:
        if redundantpool[ele['CID_scan']] == 1:
            res_new.append(ele)
    dataAll = res_new
    '''this is to distribute the weights of redundant scans'''
    for ele in dataAll:
        weightList[ele['CID_scan']] = weightList.setdefault(ele['CID_scan'], 0) + 1

    '''output current alldata'''
    # with open("{}_nonredundant.csv".format(path.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
    #     fieldnames = ['CID_scan', 'mass', 'marker', 're_rank_score', 'alpha', 'pos_a', 'a1_score', 'a2_score',
    #                   'a_mass', 'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein']
    #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #
    #     writer.writeheader()
    #     for spectrum in dataAll:
    #         a_protein = '; '.join([i for i in set(spectrum['a_protein'])])
    #         b_protein = '; '.join([i for i in set(spectrum['b_protein'])])
    #         writer.writerow(
    #             {'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'marker': spectrum['marker'],
    #              're_rank_score': spectrum['re_rank_score'], 'alpha': spectrum['alpha'],
    #              'pos_a': spectrum['pos_a'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
    #              'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
    #              'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
    #              'b_protein': b_protein, 'a2_score': spectrum['a2_score'], 'b2_score': spectrum['b2_score']})

    '''divide each type into intra and inter'''

    inter, intra = separate(dataAll)
    
    '''extra process of inter result'''
    # if inter:
    #     ave_inter_a = np.mean([ii['a1_score'] for ii in inter])
    #     ave_inter_b = np.mean([ii['b1_score'] for ii in inter])
    #     inter = [ii for ii in inter if ii['a1_score'] > 2 * ave_inter_a and ii['b1_score'] > 2 * ave_inter_b]

    res_inter = inter_fdr(inter, fdr=FDR)
    res_intra = inter_fdr(intra, fdr=FDR)

    #res_inter = sorted(res_inter, key=lambda x: -x['re_rank_score'])
    #res_intra = sorted(res_intra, key=lambda x: -x['re_rank_score'])

    '''store by two different csv files'''
    # with open("{}_intra.csv".format(path.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
    #     fieldnames = ['CID_scan', 'mass', 'marker', 'PFM_score', 'alpha', 'pos_a', 'a1_score', 'a2_score', 'a_mass',
    #                   'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein', 'q_value']
    #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #
    #     writer.writeheader()
    #     for spectrum in res_intra:
    #         a_protein = '; '.join([i for i in set(spectrum['a_protein'])])
    #         b_protein = '; '.join([i for i in set(spectrum['b_protein'])])
    #         writer.writerow(
    #             {'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'marker': spectrum['marker'],
    #              'PFM_score': spectrum['re_rank_score'], 'alpha': spectrum['alpha'],
    #              'pos_a': spectrum['pos_a'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
    #              'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
    #              'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
    #              'b_protein': b_protein, 'q_value': spectrum['q_value'],
    #              'a2_score': spectrum['a2_score'], 'b2_score': spectrum['b2_score']})
    #
    # with open("{}_inter.csv".format(path.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
    #     fieldnames = ['CID_scan', 'mass', 'marker', 'PFM_score', 'alpha', 'pos_a', 'a1_score', 'a2_score', 'a_mass',
    #                   'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein', 'q_value']
    #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #
    #     writer.writeheader()
    #     for spectrum in res_inter:
    #         a_protein = '; '.join([i.split(' ')[0] for i in set(spectrum['a_protein'])])
    #         b_protein = '; '.join([i.split(' ')[0] for i in set(spectrum['b_protein'])])
    #         writer.writerow(
    #             {'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'marker': spectrum['marker'],
    #              'PFM_score': spectrum['re_rank_score'], 'alpha': spectrum['alpha'],
    #              'pos_a': spectrum['pos_a'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
    #              'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
    #              'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
    #              'b_protein': b_protein, 'q_value': spectrum['q_value'],
    #              'a2_score': spectrum['a2_score'], 'b2_score': spectrum['b2_score']})

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
        # indivNum[ele['beta']] = indivNum.setdefault(ele['beta'], 0) + 1
        # indivNum[ele['alpha']] = indivNum.setdefault(ele['alpha'], 0) + 1
    if list(pepNum.values()):
        aveNum = np.mean(list(pepNum.values()))
    else:
        aveNum = 0
    # indivAveNum = np.mean(list(indivNum.values()))

    for ele in final_res:
        # if pepNum['{}-{}'.format(ele['alpha'], ele['beta'])] <= np.floor(aveNum / 2):
        # if pepNum[ele['alpha']] > aveNum / 2 and pepNum[ele['beta']] > aveNum / 2:
        if pepNum['{}-{}'.format(ele['alpha'], ele['beta'])] >= aveNum / 4:
            # if indivNum[ele['alpha']] > indivAveNum / 5 and indivNum[ele['beta']] > indivAveNum / 5:
            new_res.append(ele)
    final_res = sorted(new_res, key=lambda x: -x['b1_score'] - x['a1_score'])
    '''write the output file'''

    with open("{}_final.csv".format(path.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
        fieldnames = ['CID_scan', 'mass', 'marker', 'PFM_score', 'alpha', 'pos_a', 'a1_score', 'a2_score', 'a_mass',
                      'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein', 'type', 'q_value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for spectrum in final_res:
            a_protein = '; '.join([i.split(' ')[0] for i in set(spectrum['a_protein'])])
            b_protein = '; '.join([i.split(' ')[0] for i in set(spectrum['b_protein'])])
            writer.writerow(
                {'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'marker': spectrum['marker'],
                 'PFM_score': spectrum['re_rank_score'], 'alpha': spectrum['alpha'],
                 'pos_a': spectrum['pos_a'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
                 'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
                 'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
                 'b_protein': b_protein, 'q_value': spectrum['q_value'],
                 'a2_score': spectrum['a2_score'], 'b2_score': spectrum['b2_score'], 'type': spectrum['type']})

    '''this is to save inter and intra raw'''
    # with open("{}_inter_raw.csv".format(path.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
    #     fieldnames = ['CID_scan', 'mass', 'marker', 're_rank_score', 'alpha', 'pos_a', 'a1_score', 'a2_score', 'a_mass',
    #                   'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein', 'q_value']
    #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #
    #     writer.writeheader()
    #     for spectrum in inter:
    #         a_protein = spectrum['a_protein']
    #         b_protein = spectrum['b_protein']
    #         writer.writerow(
    #             {'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'marker': spectrum['marker'],
    #              're_rank_score': spectrum['re_rank_score'], 'alpha': spectrum['alpha'],
    #              'pos_a': spectrum['pos_a'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
    #              'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
    #              'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
    #              'b_protein': b_protein, 'q_value': spectrum['q_value'],
    #              'a2_score': spectrum['a2_score'], 'b2_score': spectrum['b2_score']})
    #
    # with open("{}_intra_raw.csv".format(path.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
    #     fieldnames = ['CID_scan', 'mass', 'marker', 're_rank_score', 'alpha', 'pos_a', 'a1_score', 'a2_score', 'a_mass',
    #                   'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein', 'q_value']
    #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #
    #     writer.writeheader()
    #     for spectrum in intra:
    #         a_protein = spectrum['a_protein']
    #         b_protein = spectrum['b_protein']
    #         writer.writerow(
    #             {'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'marker': spectrum['marker'],
    #              're_rank_score': spectrum['re_rank_score'], 'alpha': spectrum['alpha'],
    #              'pos_a': spectrum['pos_a'], 'a_mass': spectrum['a_mass'], 'b_mass': spectrum['b_mass'],
    #              'a1_score': spectrum['a1_score'], 'a_protein': a_protein,
    #              'beta': spectrum['beta'], 'pos_b': spectrum['pos_b'], 'b1_score': spectrum['b1_score'],
    #              'b_protein': b_protein, 'q_value': spectrum['q_value'],
    #              'a2_score': spectrum['a2_score'], 'b2_score': spectrum['b2_score']})