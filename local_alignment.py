import pickle
import numpy as np
import json
from ECL_PF import all_abs_path_of_components
import concurrent.futures


def batch_merge(ls, tol):
    """
    Merge peaks into one if they are within the tolerance.
    ls is the unordered peaks
    :param ls:
    :param tol:
    """
    ls = sorted(ls)  # ascending order
    res = [ls[0]]  # initialize peak list
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


def align(tuple_bag):
    """
    Align spectra and Construct representative spectrum to find more possible signature ions
    :param tuple_bag:
    """
    path, tol1, tol2, win = tuple_bag  # input parameters
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
    local_window = 500  # two spectra within 500 scans can be aligned
    spectra_list = all_abs_path_of_components()
    print(spectra_list)
    args = []
    for i in spectra_list:
        args.append((i, TOLERANCE1, TOLERANCE2, local_window))

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        result = executor.map(align, args)


