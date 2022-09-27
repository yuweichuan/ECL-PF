import pickle
import numpy as np
import json
from ECL_PF import all_abs_path_of_components
import concurrent.futures


def sp_similarity(l1, l2, tol):
    """ measure the similarity of CID spectrum"""
    idx1 = 0
    idx2 = 0
    count = 0
    tol = max(l1[-1][0], l2[-1][0]) * 2 * tol
    while idx1 < len(l1) and idx2 < len(l2):
        if l1[idx1][0] > l2[idx2][0] + tol:
            idx2 += 1
        elif l1[idx1][0] < l2[idx2][0] - tol:
            idx1 += 1
        else:
            count += 1
            idx1 += 1
            idx2 += 1
    return count / min(len(l1), len(l2))


def sp_merge(l1, l2, tol):
    """merge two spectrum into a big one in a sorted order"""
    if l1 is None and l2 is None:
        return None
    elif l1 is None:
        return l2
    elif l2 is None:
        return l1
    idx1 = 0
    idx2 = 0
    res = []
    tol = max(l1[-1][0], l2[-1][0]) * tol
    while idx1 < len(l1) and idx2 < len(l2):
        if l1[idx1][0] > l2[idx2][0] + tol:
            res.append(l2[idx2])
            idx2 += 1
        elif l1[idx1][0] < l2[idx2][0] - tol:
            res.append(l1[idx1])
            idx1 += 1
        else:
            mz = (l1[idx1][0] * l1[idx1][1] + l2[idx2][0] * l2[idx2][1]) / (l1[idx1][1] + l2[idx2][1])
            intense = l1[idx1][1] + l2[idx2][1]
            res.append((mz, intense))
            idx1 += 1
            idx2 += 1
    res = res + l1[idx1:] + l2[idx2:]
    return res


def batch_merge(ls, tol):
    """ls is the unordered peaks"""
    ls = sorted(ls)
    res = [ls[0]]
    start = 1
    while start < len(ls):
        if np.isclose(ls[start][0], res[-1][0], rtol=tol, atol=0):
            tempX = (ls[start][0] * ls[start][1] + res[-1][0] * res[-1][1]) / \
                    (ls[start][1] + res[-1][1])
            tempY = ls[start][1] + res[-1][1]
            res[-1] = (tempX, tempY)
        else:
            res.append(ls[start])
        start += 1
    return res


def align(tuple_bag):  # the purpose is to construct representative spectrum to find more possible signature ions
    path, tol1, tol2, win = tuple_bag
    with open(path, 'rb') as f:
        spectra = pickle.load(f)
    '''align the spectrum, group the scan number first'''
    scan1 = spectra[0]['CID_scan']
    mass1 = spectra[0]['mass']
    ls1 = spectra[0]['CID_peaks']
    group = [[[0], mass1, scan1, ls1]]  # [[scan1index, scan2,..], mass, last_scan, ls], [[scan1index, scan2,..], mass, last_scan, ls],...
    spectrumIdx = 0
    for spectrum in spectra[1:]:
        spectrumIdx += 1
        idx = -1
        flag = False
        while idx >= -len(group) and spectrum['CID_scan'] - group[idx][2] < win:
            # if np.isclose(spectrum['mass'], group[idx][1], rtol=tol1) and sp_similarity(group[idx][3], spectrum['CID_peaks'], tol2) > 0:
            if np.isclose(spectrum['mass'], group[idx][1], rtol=tol1):
                group[idx][0].append(spectrumIdx)
                group[idx][2] = spectrum['CID_scan']
                group[idx][3] = spectrum['CID_peaks']
                group.append(group.pop(idx))
                flag = True
                break
            else:
                idx -= 1
        if not flag:
            curScan = spectrum['CID_scan']
            curMass = spectrum['mass']
            curLs = spectrum['CID_peaks']
            group.append([[spectrumIdx], curMass, curScan, curLs])
    res = []
    pc = 0
    for sub in group:
        pc += 1
        cpeak = []
        for scanIdx in sub[0]:
            cpeak += spectra[scanIdx]['CID_peaks']

        cpeak = batch_merge(cpeak, tol2)

        for scanIdx in sub[0]:
            res.append({'mass': spectra[scanIdx]['mass'], 'CID_scan': spectra[scanIdx]['CID_scan'],
                        'charge': spectra[scanIdx]['charge'], 'rep_CID_peaks': cpeak, 'repNum': len(sub[0]),
                        'CID_peaks': spectra[scanIdx]['CID_peaks'], 'ETD_peaks': spectra[scanIdx]['ETD_peaks']})
    print('Done! {}'.format(path))
    with open('{}'.format(path), 'wb') as f:
        pickle.dump(res, f)


if __name__ == '__main__':
    with open('ECLPF_conf', 'r') as conf_file:
        conf = json.load(conf_file)
    TOLERANCE1 = conf['ms1_tol']
    TOLERANCE2 = conf['ms2_tol']
    local_window = 500
    spectra_list = all_abs_path_of_components()
    print(spectra_list)
    args = []
    for i in spectra_list:
        args.append((i, TOLERANCE1, TOLERANCE2, local_window))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        result = executor.map(align, args)


