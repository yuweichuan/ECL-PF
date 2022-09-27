from fragment import fast_b_y_fragment, fast_b_y_c_z_fragment
from precursor_discovery import find_precursor_ab_mass, extract_peptide
from pyteomics import parser, mass
import numpy as np
import pickle
import os
import time
import copy
import math


def spec_pre_process(mz_list):  # sorted [(mass1,intensity1),(mass2,intensity2),(),()...]
    """spectrum normalization in local range (100 Da) since the weighted average amino acid mass is 111.11 Da,
    return sorted [(mass1,norm_intensity1),(mass2,norm_intensity2),(),()...]"""
    mz_int_list = list(mz_list)
    min_mz = mz_int_list[0][0]
    max_mz = mz_int_list[-1][0]
    res = []
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


def cid_x_corr_score(mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass, xl_site,
                     proN=True):  # tol ppm
    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:
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
        xcorr = 0
        theo_list = fast_b_y_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)
        index_theo = 0
        index_exp = 0
        while index_exp < len(mz_int_list) and index_theo < len(theo_list):
            if mz_int_list[index_exp][0] - theo_list[index_theo] < \
                    - tol * max(mz_int_list[index_exp][0], theo_list[index_theo]):
                index_exp += 1
            elif mz_int_list[index_exp][0] - theo_list[index_theo] > \
                    tol * max(mz_int_list[index_exp][0], theo_list[index_theo]):
                index_theo += 1
            else:
                xcorr += mz_int_list[index_exp][1]
                index_exp += 1
                index_theo += 1
        res.append((xcorr, peptide_sequence, pos_ele))
    return max(res)


def cid_etd_x_corr_score(cid_mz_int_list, etd_mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass,
                         xl_site, proN=True):  # tol ppm
    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:
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
        cid_theo_list, etd_theo_list = \
            fast_b_y_c_z_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)
        cid_xcorr = 0
        index_theo = 0
        index_exp = 0
        while index_exp < len(cid_mz_int_list) and index_theo < len(cid_theo_list):
            if cid_mz_int_list[index_exp][0] - cid_theo_list[index_theo] < \
                    - tol * max(cid_mz_int_list[index_exp][0], cid_theo_list[index_theo]):
                index_exp += 1
            elif cid_mz_int_list[index_exp][0] - cid_theo_list[index_theo] > \
                    tol * max(cid_mz_int_list[index_exp][0], cid_theo_list[index_theo]):
                index_theo += 1
            else:
                cid_xcorr += cid_mz_int_list[index_exp][1]
                index_exp += 1
                index_theo += 1

        etd_xcorr = 0
        index_theo = 0
        index_exp = 0
        while index_exp < len(etd_mz_int_list) and index_theo < len(etd_theo_list):
            if etd_mz_int_list[index_exp][0] - etd_theo_list[index_theo] < \
                    - tol * max(etd_mz_int_list[index_exp][0], etd_theo_list[index_theo]):
                index_exp += 1
            elif etd_mz_int_list[index_exp][0] - etd_theo_list[index_theo] > \
                    tol * max(etd_mz_int_list[index_exp][0], etd_theo_list[index_theo]):
                index_theo += 1
            else:
                etd_xcorr += etd_mz_int_list[index_exp][1]
                index_exp += 1
                index_theo += 1

        res.append(((cid_xcorr + etd_xcorr) / 2, peptide_sequence, pos_ele))
    return max(res)


def cid_p_score(mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass, xl_site, xlMass,
                proN=True):  # tol ppm
    """return -log10(p_score)"""
    x = 1 / 111.1 * 2 * (pre_mass * tol * 2)
    factor = chain_mass / (pre_mass - xlMass) * len(mz_int_list)
    # replace sequence length with mass
    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:
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
        theo_list = fast_b_y_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)
        index_theo = 0
        index_exp = 0
        while index_exp < len(mz_int_list) and index_theo < len(theo_list):
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
    """return -log10(p_score)"""
    x = 1 / 111.1 * 4 * (pre_mass * tol * 2)
    factor = chain_mass / (pre_mass - xlMass) * (len(cid_mz_int_list) + len(etd_mz_int_list))
    # replace sequence length with mass

    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:
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
            fast_b_y_c_z_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)
        index_theo = 0
        index_exp = 0
        while index_exp < len(cid_mz_int_list) and index_theo < len(cid_theo_list):
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
        while index_exp < len(etd_mz_int_list) and index_theo < len(etd_theo_list):
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

        n = int(n / 2)
        for ii in range(n):
            p_score -= math.exp(-x * factor) * (x * factor) ** ii / math.factorial(ii)
        if p_score < 0:
            p_score = 1E-20

        res.append((-math.log10(p_score), peptide_sequence, pos_ele))
    return max(res)


def cid_merox_score(mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass, xl_site, ml, ms, charge,
                    proN=True):
    """input raw spectrum, modified merox score without p index and ignore Fxl because it is always zero"""
    d = len(mz_int_list)  # number of data points
    if d > 460:
        d = 460
    if d < 140:
        d = 140
    maxIntensity = max([ii[1] for ii in mz_int_list])
    fTotal = sum([ii[1] for ii in mz_int_list])
    cutoff = 0.1 * maxIntensity
    abSig = len([ii[1] for ii in mz_int_list if ii[1] > cutoff])  # number of abundant signals (above cutoff)

    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:
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

    '''calculate the matched ions'''
    for pos_ele in pos:
        h = 0  # number of identified ions
        Ik = 0  # abundance of identified signals
        k = 0  # number of identified abundant signals
        theo_list = fast_b_y_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)
        index_theo = 0
        index_exp = 0
        while index_exp < len(mz_int_list) and index_theo < len(theo_list):
            if mz_int_list[index_exp][0] - theo_list[index_theo] < \
                    - tol * max(mz_int_list[index_exp][0], theo_list[index_theo]):
                index_exp += 1
            elif mz_int_list[index_exp][0] - theo_list[index_theo] > \
                    tol * max(mz_int_list[index_exp][0], theo_list[index_theo]):
                index_theo += 1
            else:
                h += 1
                Ik += mz_int_list[index_exp][1]
                k += 1 if mz_int_list[index_exp][1] > cutoff else 0
                index_exp += 1
                index_theo += 1
        fInt = Ik / fTotal  # factor for signal intensity
        score = -50 * math.log10(0.1 * np.exp(-k / abSig) + 0.05 * np.exp(-20 * h / d) + 0.1 * np.exp(-k / 6) +
                                 0.3 / charge * 20 ** (-fInt ** 2)) + 50 * math.log10(0.25 + 0.3 / charge)

        res.append((score, peptide_sequence, pos_ele))
    return max(res)


def cid_etd_merox_score(cid_mz_int_list, etd_mz_int_list, peptide_sequence, pre_mass, chain_mass, tol, fra_aa_mass,
                        xl_site, ml, ms, charge, proN=True):  # tol ppm
    """input raw spectrum, modified merox score without p index and do not distinguish charge status in Fxl"""
    d = len(cid_mz_int_list) + len(etd_mz_int_list)  # number of data points
    if d > 460:
        d = 460
    if d < 140:
        d = 140
    cidMaxIntensity = max([ii[1] for ii in cid_mz_int_list])
    etdMaxIntensity = max([ii[1] for ii in etd_mz_int_list])
    cid_cutoff = 0.1 * cidMaxIntensity
    etd_cutoff = 0.1 * etdMaxIntensity
    fTotal = sum([ii[1] for ii in cid_mz_int_list]) + sum([ii[1] for ii in etd_mz_int_list])
    abSig = len([ii[1] for ii in cid_mz_int_list if ii[1] > 0.1 * cidMaxIntensity]) \
        + len([ii[1] for ii in etd_mz_int_list if ii[1] > 0.1 * etdMaxIntensity])

    proteinN = False
    proteinC = False
    if '[' == peptide_sequence[0]:
        peptide_sequence = peptide_sequence[1:]
        proteinN = proN
    if ']' == peptide_sequence[-1]:
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
        h = 0  # number of identified ions
        Ik = 0  # abundance of identified signals
        k = 0  # number of identified abundant signals
        cid_theo_list, etd_theo_list = \
            fast_b_y_c_z_fragment(peptide_sequence, pos_ele, chain_mass, pre_mass, fra_aa_mass)
        index_theo = 0
        index_exp = 0
        while index_exp < len(cid_mz_int_list) and index_theo < len(cid_theo_list):
            if cid_mz_int_list[index_exp][0] - cid_theo_list[index_theo] < \
                    - tol * max(cid_mz_int_list[index_exp][0], cid_theo_list[index_theo]):
                index_exp += 1
            elif cid_mz_int_list[index_exp][0] - cid_theo_list[index_theo] > \
                    tol * max(cid_mz_int_list[index_exp][0], cid_theo_list[index_theo]):
                index_theo += 1
            else:
                h += 1
                Ik += cid_mz_int_list[index_exp][1]
                k += 1 if cid_mz_int_list[index_exp][1] > cid_cutoff else 0
                index_exp += 1
                index_theo += 1

        index_theo = 0
        index_exp = 0
        while index_exp < len(etd_mz_int_list) and index_theo < len(etd_theo_list):
            if etd_mz_int_list[index_exp][0] - etd_theo_list[index_theo] < \
                    - tol * max(etd_mz_int_list[index_exp][0], etd_theo_list[index_theo]):
                index_exp += 1
            elif etd_mz_int_list[index_exp][0] - etd_theo_list[index_theo] > \
                    tol * max(etd_mz_int_list[index_exp][0], etd_theo_list[index_theo]):
                index_theo += 1
            else:
                h += 1
                Ik += etd_mz_int_list[index_exp][1]
                k += 1 if etd_mz_int_list[index_exp][1] > etd_cutoff else 0
                index_exp += 1
                index_theo += 1

        fInt = Ik / fTotal
        score = -50 * math.log10(0.1 * np.exp(-k / abSig) + 0.05 * np.exp(-20 * h / d) + 0.1 * np.exp(-k / 6) +
                                 0.3 / charge * 20 ** (-fInt ** 2)) + 50 * math.log10(0.25 + 0.3 / charge)
        res.append((score, peptide_sequence, pos_ele))
    return max(res)


if __name__ == '__main__':
    Link_site = ['K']  # support multiple sites such as ['K', 'R']
    Fix_mod = {'car': [57.021464, ['C']]}
    Var_mod = {'28': [28.0313, ['K', 'Peptide-nterm']],
               '34': [34.0631, ['K', 'Peptide-nterm']],
               'oxi': [15.99, ['M']]}

    '''add modification masses into the dict'''
    new_mass = {}
    for i, j in Fix_mod.items():
        new_mass[i] = j[0]

    for i, j in Var_mod.items():
        new_mass[i] = j[0]

    AA_mass = dict(mass.std_aa_mass)
    AA_mass.update(new_mass)
    xl_mass = 509.0963
    large_mass = 455.0868
    small_mass = 54.0106
    TOL = 2e-5
    pre_TOL = 3e-5
    with open(r"D:\OneDrive - HKUST Connect\Research\ECL_X\ECLX_src\data\pickled_pair_wise\MS180768_L1_hcd26_etd",
              'rb') as f:
        a = pickle.load(f)
    for i in a:
        if i['CID_scan'] == 8898:
            mz_int = [ii for ii in i['CID_peaks'] if ii[0] < i['mass']]
            print(mz_int)
            PRECURSOR_MASS = i['mass']
    TT = 1
    start1 = time.time()
    for tic in range(TT):
        ab = find_precursor_ab_mass(mz_int, PRECURSOR_MASS, xl_mass, large_mass, small_mass, TOL)
        if ab[0][2] > 2:
            ab = [ele for ele in ab if ele[2] > 2]
    print(mz_int)
    print(len(ab))
    start2 = time.time()
    for tic in range(TT):
        processed_mz_int = spec_pre_process(copy.copy(mz_int))
        Xcorr = 0
    start3 = time.time()
    sub_ex = 0
    sub_sc = 0
    for tic in range(TT):
        for ab_ele in ab:
            tic1 = time.time()
            a_xcorr = 0
            b_xcorr = 0
            alpha = ab_ele[1]
            beta = ab_ele[0]
            alpha_candidate = extract_peptide(alpha, alpha * pre_TOL)
            beta_candidate = extract_peptide(beta, beta * pre_TOL)
            toc1 = time.time()
            sub_ex += toc1 - tic1
            tic2 = time.time()
            if alpha_candidate and beta_candidate:
                for alpha_ele in alpha_candidate:
                    alpha_res = \
                        cid_x_corr_score(processed_mz_int, alpha_ele[1], PRECURSOR_MASS, alpha_ele[0], TOL, AA_mass,
                                         Link_site)
                    if alpha_res[0] >= a_xcorr:
                        pepa = alpha_res[1]
                        a_xcorr = alpha_res[0]
                for beta_ele in beta_candidate:
                    beta_res = \
                        cid_x_corr_score(processed_mz_int, beta_ele[1], PRECURSOR_MASS, beta_ele[0], TOL, AA_mass,
                                         Link_site)
                    if beta_res[0] >= b_xcorr:
                        pepb = beta_res[1]
                        b_xcorr = beta_res[0]
            if b_xcorr + a_xcorr > Xcorr:
                pep = [pepb, pepa]
                Xcorr = b_xcorr + a_xcorr
            toc2 = time.time()
            sub_sc += toc2 - tic2
    print(9.540314435958862, 0.9406578540802002, 5.416500091552734, 8.658490419387817, 0.8662509918212891)
    print(time.time() - start3, start3 - start2, start2 - start1, sub_ex, sub_sc)
    try:
        print(pep)
    except:
        print('no pep')
    strat = time.time()
