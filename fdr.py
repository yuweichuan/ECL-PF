import json
import csv
import numpy as np


def flatten(ele):
    if not isinstance(ele, (list, )):
        return [ele]
    else:
        b = []
        for item in ele:
            b += flatten(item)
    return b


def fdr(list_dict):
    """considering ETD and CID"""
    """"[{'CID_scan':_, 'mass':_, 'alpha1':[score, seq, pos, des, mass, ETD score],
    'beta1':[score, seq, pos, des, mass, ETD score],
    'alpha2':score,
    'beta2':score,
    'marker':_,},{},{}]
    return sorted list dict by q value"""
    for ele in list_dict:  # delete the decoy if target exist

        ele['alpha1'][3] = ele['alpha1'][3].split('$')
        ele['alpha1'][3] = [int(alpha1_ele) for alpha1_ele in ele['alpha1'][3]]
        ele['beta1'][3] = ele['beta1'][3].split('$')
        ele['beta1'][3] = [int(beta1_ele) for beta1_ele in ele['beta1'][3]]

        TF = []
        for des in ele['alpha1'][3]:
            TF.append(not bool(des % 2))
        if any(TF) and not all(TF):
            new = list(ele['alpha1'][3])
            for tf in range(len(TF)):
                if TF[tf]:
                    new.remove(ele['alpha1'][3][tf])
            ele['alpha1'][3] = new
        TF = []
        for des in ele['beta1'][3]:
            TF.append(not bool(des % 2))
        if any(TF) and not all(TF):
            new = list(ele['beta1'][3])
            for tf in range(len(TF)):
                if TF[tf]:
                    new.remove(ele['beta1'][3][tf])
            ele['beta1'][3] = new
    list_dict = [ele for ele in list_dict if ele['alpha1'][0] + ele['beta1'][0] > 0]
    # delete negative and zero samples, therefore remove sequence 'PEPTIDE'
    '''define the re-rank score first'''
    for ele in list_dict:
        ele['re_rank'] = ele['protein_score']

    list_dict = \
        sorted(list_dict, key=lambda x: -x['re_rank'])
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict) - 1, -1, -1):
        if list_dict[idx]['fdr'] < minfdr:
            minfdr = list_dict[idx]['fdr']
        list_dict[idx]['q_value'] = minfdr
    return list_dict


def specific_separate_fdr(list_dict, q_val):
    """considering ETD and CID, choose the greater one"""
    """"[{'CID_scan':_, 'mass':_, 'alpha1':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'beta1':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'alpha2':score,
    'beta2':score,
    'marker':_,},{},{}]
    return sorted list dict by q value"""
    for ele in list_dict:  # delete the decoy if target exist
        TF = []
        for des in ele['alpha1'][3]:
            TF.append(not bool(des % 2))
        if any(TF) and not all(TF):
            new = list(ele['alpha1'][3])
            for tf in range(len(TF)):
                if TF[tf]:
                    new.remove(ele['alpha1'][3][tf])
            ele['alpha1'][3] = new
        TF = []
        for des in ele['beta1'][3]:
            TF.append(not bool(des % 2))
        if any(TF) and not all(TF):
            new = list(ele['beta1'][3])
            for tf in range(len(TF)):
                if TF[tf]:
                    new.remove(ele['beta1'][3][tf])
            ele['beta1'][3] = new
    list_dict = [ele for ele in list_dict if ele['alpha1'][0] > 0 and ele['beta1'][0] > 0]
    # delete negative and zero samples
    '''group the data into different marker'''
    list_dict4 = [ele for ele in list_dict if ele['marker'] == 4]
    list_dict3 = [ele for ele in list_dict if ele['marker'] == 3]
    list_dict2 = [ele for ele in list_dict if ele['marker'] == 2]
    list_dict15 = [ele for ele in list_dict if ele['marker'] == 1.5]
    list_dict125 = [ele for ele in list_dict if ele['marker'] == 1.25]
    list_dict1 = [ele for ele in list_dict if ele['marker'] == 1]
    '''process list_dict4'''
    list_dict4 = \
        sorted(list_dict4, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict4:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict4) - 1, -1, -1):
        if list_dict4[idx]['fdr'] < minfdr:
            minfdr = list_dict4[idx]['fdr']
        list_dict4[idx]['q_value'] = minfdr
        list_dict4[idx]['re_rank'] = \
            max(list_dict4[idx]['alpha1'][0], list_dict4[idx]['alpha1'][5]) * \
            max(list_dict4[idx]['beta1'][0], list_dict4[idx]['beta1'][5])
    list_dict4 = [ele for ele in list_dict4 if ele['q_value'] <= q_val and bool(ele['alpha1'][3][0] % 2) and
                  bool(ele['beta1'][3][0] % 2)]

    '''process list_dict3'''
    list_dict3 = \
        sorted(list_dict3, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict3:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict3) - 1, -1, -1):
        if list_dict3[idx]['fdr'] < minfdr:
            minfdr = list_dict3[idx]['fdr']
        list_dict3[idx]['q_value'] = minfdr
        list_dict3[idx]['re_rank'] = \
            max(list_dict3[idx]['alpha1'][0], list_dict3[idx]['alpha1'][5]) * \
            max(list_dict3[idx]['beta1'][0], list_dict3[idx]['beta1'][5])
    list_dict3 = [ele for ele in list_dict3 if ele['q_value'] <= q_val and bool(ele['alpha1'][3][0] % 2) and
                  bool(ele['beta1'][3][0] % 2)]

    '''process list_dict2'''
    list_dict2 = \
        sorted(list_dict2, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict2:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict2) - 1, -1, -1):
        if list_dict2[idx]['fdr'] < minfdr:
            minfdr = list_dict2[idx]['fdr']
        list_dict2[idx]['q_value'] = minfdr
        list_dict2[idx]['re_rank'] = \
            max(list_dict2[idx]['alpha1'][0], list_dict2[idx]['alpha1'][5]) * \
            max(list_dict2[idx]['beta1'][0], list_dict2[idx]['beta1'][5])
    list_dict2 = [ele for ele in list_dict2 if ele['q_value'] <= q_val and bool(ele['alpha1'][3][0] % 2) and
                  bool(ele['beta1'][3][0] % 2)]

    '''process list_dict15'''
    list_dict15 = \
        sorted(list_dict15, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict15:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict15) - 1, -1, -1):
        if list_dict15[idx]['fdr'] < minfdr:
            minfdr = list_dict15[idx]['fdr']
        list_dict15[idx]['q_value'] = minfdr
        list_dict15[idx]['re_rank'] = \
            max(list_dict15[idx]['alpha1'][0], list_dict15[idx]['alpha1'][5]) * \
            max(list_dict15[idx]['beta1'][0], list_dict15[idx]['beta1'][5])
    list_dict15 = [ele for ele in list_dict15 if ele['q_value'] <= q_val and bool(ele['alpha1'][3][0] % 2) and
                   bool(ele['beta1'][3][0] % 2)]
    '''process list_dict125'''
    list_dict125 = \
        sorted(list_dict125, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict125:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict125) - 1, -1, -1):
        if list_dict125[idx]['fdr'] < minfdr:
            minfdr = list_dict125[idx]['fdr']
        list_dict125[idx]['q_value'] = minfdr
        list_dict125[idx]['re_rank'] = \
            max(list_dict125[idx]['alpha1'][0], list_dict125[idx]['alpha1'][5]) * \
            max(list_dict125[idx]['beta1'][0], list_dict125[idx]['beta1'][5])
    list_dict125 = [ele for ele in list_dict125 if ele['q_value'] <= q_val and bool(ele['alpha1'][3][0] % 2) and
                    bool(ele['beta1'][3][0] % 2)]

    '''process list_dict1'''
    list_dict1 = \
        sorted(list_dict1, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    TT = 0
    TD = 0
    DD = 0
    for ele in list_dict1:
        tga = not bool(ele['alpha1'][3][0] % 2)
        tgb = not bool(ele['beta1'][3][0] % 2)
        if tgb and tga:
            DD += 1
        elif not tga and not tgb:
            TT += 1
        else:
            TD += 1
        if TT == 0:
            ele['fdr'] = 1
        else:
            if TD > DD:
                ele['fdr'] = (TD - DD) / TT
            else:
                ele['fdr'] = DD / TT
    '''trace back'''
    minfdr = 1
    for idx in range(len(list_dict1) - 1, -1, -1):
        if list_dict1[idx]['fdr'] < minfdr:
            minfdr = list_dict1[idx]['fdr']
        list_dict1[idx]['q_value'] = minfdr
        list_dict1[idx]['re_rank'] = \
            max(list_dict1[idx]['alpha1'][0], list_dict1[idx]['alpha1'][5]) * \
            max(list_dict1[idx]['beta1'][0], list_dict1[idx]['beta1'][5])
    list_dict1 = [ele for ele in list_dict1 if ele['q_value'] <= q_val and bool(ele['alpha1'][3][0] % 2) and
                  bool(ele['beta1'][3][0] % 2)]
    res = list_dict4 + list_dict3 + list_dict2 + list_dict15 + list_dict125 + list_dict1
    res = sorted(res, key=lambda x: (-max(x['alpha1'][0], x['alpha1'][5]) * max(x['beta1'][0], x['beta1'][5])))
    return res


if __name__ == '__main__':
    with open('psm', 'r') as f:
        spectra = json.load(f)
    a = fdr(spectra)
    with open(f"{'ECLX'}.csv", 'w', newline='') as csvfile:
        fieldnames = ['CID_scan', 'mass', 'alpha', 'pos_a', 'a_score', 'beta', 'pos_b', 'b_score', 'q_value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for spectrum in a:
            writer.writerow({'CID_scan': spectrum['CID_scan'], 'mass': spectrum['mass'], 'alpha': spectrum['alpha1'][1],
                             'pos_a': spectrum['alpha1'][2], 'a_score': spectrum['alpha1'][0],
                             'beta': spectrum['beta1'][1], 'pos_b': spectrum['beta1'][2],
                             'b_score': spectrum['beta1'][0], 'q_value': spectrum['q_value']})
