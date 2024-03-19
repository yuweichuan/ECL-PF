import sys
import csv
import numpy as np
import os
import json


def separate(data):
    """
    Separate the ran out data as inter and intra-protein groups.
    Inter group includes homo-oligomers
    :param data:
    """
    inter = []  # inter-protein group
    intra = []  # intra-protein group

    for sub_data in data:
        inter.append(dict(sub_data))
        if set(sub_data['a_protein']) & set(sub_data['b_protein']):
            intra.append(dict(sub_data))

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

    '''Trace back the q-value'''
    minfdr = 1
    for idx in range(len(data) - 1, -1, -1):
        if data[idx]['fdr'] < minfdr:
            minfdr = data[idx]['fdr']
        data[idx]['q_value'] = minfdr

    data = [sub_data for sub_data in data if sub_data['q_value'] <= fdr and 'DECOY' not in sub_data['a_protein'][0] and
            'DECOY' not in sub_data['b_protein'][0]]
    return data


if __name__ == '__main__':
    """The main function to separate the data into intra 
    and inter-protein groups and calculate the fdr"""
    with open('ECLPF_conf', 'r') as conf_file:  # parameters file
        conf = json.load(conf_file)

    dataAll = []
    path = sys.argv[1]
    outpath = sys.argv[2]
    FDR = float(sys.argv[3])
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

    'Sort CSMs by scan number'
    dataAll = sorted(dataAll, key=lambda x: (-x['CID_scan'], -x['b1_score'] - x['a1_score']))

    '''Only retain the highest scored ones for each scan number'''
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

    '''Remove redundant decoys if the target candidates exists for the same peptide'''
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

    '''Distinguish ambiguity by related ones if we can.'''
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

    '''Remove redundant (optional)'''
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

    '''Divide each type into intra and inter'''
    inter, intra = separate(dataAll)
    
    '''Calculate FDR separately'''
    res_inter = inter_fdr(inter, fdr=FDR)
    res_intra = inter_fdr(intra, fdr=FDR)

    '''Post processing'''
    final_res = []
    inter_scan = dict()
    res_index = 0

    for ele in res_inter:
        ele['type'] = 'inter'
        final_res.append(ele)
        inter_scan[ele['CID_scan']] = res_index
        res_index += 1
    for ele in res_intra:  # if FDR conflict in intra and inter group, take the min one
        if ele['CID_scan'] in inter_scan.keys():
            if len(final_res[inter_scan[ele['CID_scan']]]['a_protein']) == 1 and\
                    final_res[inter_scan[ele['CID_scan']]]['a_protein'] == final_res[inter_scan[ele['CID_scan']]]['b_protein']:

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

    '''Calculate average number of one cross-linked peptides being identified and remove less confident'''
    new_res = []
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

    '''Write the output file'''
    with open(outpath, 'w', newline='') as csvfile:
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

    os.remove(path)
