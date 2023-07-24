from fragment import fast_b_y_fragment, fast_b_y_c_z_fragment
import numpy as np
import math


def spec_pre_process(mz_list):  # sorted [(mass1,intensity1),(mass2,intensity2),(),()...]
    """
    Spectrum normalization in SEQUEST pre-processing way
    :param mz_list:
    """
    mz_int_list = list(mz_list)
    min_mz = mz_int_list[0][0]  # min m/z in the spectrum
    max_mz = mz_int_list[-1][0]  # max m/z in the spectrum
    res = []  # returned the peak list
    for start in range(int((min_mz // 100) * 100), int(np.ceil(max_mz + 0.0001)),
                       100):  # 0.0001 leeway for extreme case
        local_max = 0
        local_sub = []
        end = start + 100
        while mz_int_list and mz_int_list[0][0] < end:
            local_sub.append(mz_int_list.pop(0))
        for ele in local_sub:
            local_max = max(local_max, ele[1])
        if local_max == 0:
            continue
        else:
            normalized_factor = 50 / local_max ** 0.5
            for ele in local_sub:
                res.append((ele[0], round(ele[1] ** 0.5 * normalized_factor, 3)))
    return res


def cid_p_score(mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass, xl_site, xlMass,
                proN=True):  # tolerance in ppm unit
    """
    XlinkX scoring function to match the spectrum with peptides.
    return -log10(p_score)
    :param mz_int_list:
    :param peptide_sequence:
    :param pre_mass:
    :param chain_mass:
    :param tol:
    :param fra_aa_mass:
    :param xl_site:
    :param xlMass:
    :param proN:
    """
    x = 1 / 111.1 * 2 * (pre_mass * tol * 2)
    factor = chain_mass / (pre_mass - xlMass) * len(mz_int_list)
    # replace sequence length with mass
    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:  # consider n-term
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:  # consider c-term
        peptide_sequence = peptide_sequence[:-1]
        proteinC = True

    '''Start to combine the peptides and generate theoretical peaks'''
    mid_pos = [mid_ele for mid_ele in range(len(peptide_sequence)) if str.isupper(peptide_sequence[mid_ele])]
    markList = [0] * (len(mid_pos) - 1)
    if proteinN:
        markList[0] += 1

    if peptide_sequence[mid_pos[0]] in xl_site:
        markList[0] += 2

    markList[0] = 2 if markList[0] > 2 else markList[0]

    if mid_pos[0] > 0:
        if ')(' in peptide_sequence[:mid_pos[0]]:
            markList[0] -= 2
        else:
            markList[0] -= 1

    for mark_ele in range(1, len(markList)):
        if peptide_sequence[mid_pos[mark_ele]] in xl_site:
            markList[mark_ele] += 1
        if mid_pos[mark_ele] - mid_pos[mark_ele - 1] > 1:
            markList[mark_ele] -= 1
    pos = [idx for idx in range(len(markList)) if markList[idx] > 0]
    if peptide_sequence[-1] in xl_site and proteinC:
        pos.append(len(markList))

    '''no link site in the sequence'''
    if not pos:
        return None
    res = []
    for pos_ele in pos:  # possible positions of the cross-linked sites
        p_score = 1  # initialize the p-score
        n = 0
        theo_list = fast_b_y_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)  # generate theoretical peaks
        index_theo = 0
        index_exp = 0
        while index_exp < len(mz_int_list) and index_theo < len(theo_list):  # match the peaks
            if mz_int_list[index_exp][0] - theo_list[index_theo] < \
                    - tol * max(mz_int_list[index_exp][0], theo_list[index_theo]):
                index_exp += 1
            elif mz_int_list[index_exp][0] - theo_list[index_theo] > \
                    tol * max(mz_int_list[index_exp][0], theo_list[index_theo]):
                index_theo += 1
            else:
                p_score -= math.exp(-x * factor) * (x * factor) ** n / math.factorial(n)
                n += 1
                index_exp += 1
                index_theo += 1
        if p_score < 0:
            p_score = 1E-20
        res.append((-math.log10(p_score), peptide_sequence, pos_ele))
    return max(res)


def cid_etd_p_score(cid_mz_int_list, etd_mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass,
                    xl_site, xlMass, proN=True):  # tol ppm
    """
    XlinkX scoring function to match the spectrum with peptides. b/y c/z ions
    return -log10(p_score)
    :param cid_mz_int_list:
    :param etd_mz_int_list:
    :param peptide_sequence:
    :param pre_mass:
    :param chain_mass:
    :param tol:
    :param fra_aa_mass:
    :param xl_site:
    :param xlMass:
    :param proN:
    """
    x = 1 / 111.1 * 4 * (pre_mass * tol * 2)
    factor = chain_mass / (pre_mass - xlMass) * (len(cid_mz_int_list) + len(etd_mz_int_list))
    # replace sequence length with mass
    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:  # consider n-term
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:  # consider c-term
        peptide_sequence = peptide_sequence[:-1]
        proteinC = True

    mid_pos = [mid_ele for mid_ele in range(len(peptide_sequence)) if str.isupper(peptide_sequence[mid_ele])]
    markList = [0] * (len(mid_pos) - 1)
    if proteinN:
        markList[0] += 1

    if peptide_sequence[mid_pos[0]] in xl_site:
        markList[0] += 2

    markList[0] = 2 if markList[0] > 2 else markList[0]

    if mid_pos[0] > 0:
        if ')(' in peptide_sequence[:mid_pos[0]]:
            markList[0] -= 2
        else:
            markList[0] -= 1

    for mark_ele in range(1, len(markList)):
        if peptide_sequence[mid_pos[mark_ele]] in xl_site:
            markList[mark_ele] += 1
        if mid_pos[mark_ele] - mid_pos[mark_ele - 1] > 1:
            markList[mark_ele] -= 1
    pos = [idx for idx in range(len(markList)) if markList[idx] > 0]
    if peptide_sequence[-1] in xl_site and proteinC:
        pos.append(len(markList))

    '''no link site in the sequence'''
    if not pos:
        return None
    res = []
    for pos_ele in pos:
        p_score = 1
        n = 0
        cid_theo_list, etd_theo_list = \
            fast_b_y_c_z_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)  # theoretical b/y c/z peaks
        index_theo = 0
        index_exp = 0
        while index_exp < len(cid_mz_int_list) and index_theo < len(cid_theo_list):  # match b/y peaks
            if cid_mz_int_list[index_exp][0] - cid_theo_list[index_theo] < \
                    - tol * max(cid_mz_int_list[index_exp][0], cid_theo_list[index_theo]):
                index_exp += 1
            elif cid_mz_int_list[index_exp][0] - cid_theo_list[index_theo] > \
                    tol * max(cid_mz_int_list[index_exp][0], cid_theo_list[index_theo]):
                index_theo += 1
            else:
                n += 1
                index_exp += 1
                index_theo += 1

        index_theo = 0
        index_exp = 0
        while index_exp < len(etd_mz_int_list) and index_theo < len(etd_theo_list):  # match c/z peaks
            if etd_mz_int_list[index_exp][0] - etd_theo_list[index_theo] < \
                    - tol * max(etd_mz_int_list[index_exp][0], etd_theo_list[index_theo]):
                index_exp += 1
            elif etd_mz_int_list[index_exp][0] - etd_theo_list[index_theo] > \
                    tol * max(etd_mz_int_list[index_exp][0], etd_theo_list[index_theo]):
                index_theo += 1
            else:
                n += 1
                index_exp += 1
                index_theo += 1

        n = int(n / 2)  # divided matched case by 2 to make fair comparison to CID alone case
        for ii in range(n):
            p_score -= math.exp(-x * factor) * (x * factor) ** ii / math.factorial(ii)
        if p_score < 0:
            p_score = 1E-20

        res.append((-math.log10(p_score), peptide_sequence, pos_ele))
    return max(res)

