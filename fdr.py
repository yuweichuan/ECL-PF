def fdr(list_dict):
    """
    A coarse fdr that will deprecate in the final step. Will be removed in the future.
    [{'CID_scan':_, 'mass':_, 'alpha1':[score, seq, pos, des, mass, ETD score],
    'beta1':[score, seq, pos, des, mass, ETD score],
    'alpha2':score,
    'beta2':score,
    'marker':_,},{},{}]
    return sorted list dict by q value
    :param list_dict:
    """
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
