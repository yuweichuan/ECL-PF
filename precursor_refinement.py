from pyteomics import mzxml
import numpy as np
from scipy.stats import binom
import pickle
from collections import deque
import concurrent.futures
import json
import os
from bisect import bisect_left, bisect_right


def pick_peak(lis, preMz, windowSize):
    """
    Pick peaks inside the window
    :param lis:
    :param preMz:
    :param windowSize:
    """
    lower = bisect_left(lis[0], preMz - windowSize / 2)  # window left bound
    upper = bisect_right(lis[0], preMz + windowSize / 2)  # window right bound
    res = list(zip(lis[0][lower:upper], lis[1][lower:upper]))
    res = [peak for peak in res if peak[1] > 0]  # remove zero samples if exists
    return res


def binomial(_mass):
    """
    Binomial distribution of the theoretical isotope cluster.
    Return the tuple of the first nth element until the last
    probability less than 1e-2 left
    :param _mass:
    """
    prob = 0.00062567  # averagine model, the empirical probability of the occurrence of a isotope mass
    res = []  # output theoretical peak intensities
    accumulation = 0  # intensity accumulation
    iso = 0  # isotope start position
    while accumulation < 1 - 1e-2:
        ele = binom.pmf(iso, int(_mass), prob)
        accumulation += ele
        res.append(ele)
        iso += 1
    max_ratio = 0
    for index in range(1, len(res)):
        if max(res[index - 1], res[index]) / min(res[index - 1], res[index]) > max_ratio:
            max_ratio = max(res[index - 1], res[index]) / min(res[index - 1], res[index])

    return res, max_ratio  # max ratio, in case of the noise or other overlapped cluster


def top_res_merge(cls, fls, ms1tol):  # current res lis, final res lis
    """
    Merge current calibrated mass with final calibrated mass list, delete the same ones
    :param cls:
    :param fls:
    :param ms1tol:
    """
    cur = cls + fls  # current result
    cur.sort(key=lambda x: x[1])  # sort by intensity
    idx = 1
    res = [cur[0]]
    while idx < len(cur):
        if np.isclose(cur[idx][1], res[-1][1], rtol=2 * ms1tol, atol=0):
            res[-1] = (max(res[-1][0], cur[idx][0]), res[-1][1])
        else:
            res.append(cur[idx])
        idx += 1
    res = sorted(res, reverse=True)
    return res


def precursor_determination(base_lis, ref_lis, exp_lis, preCharge, ms1tol):  # exp_lis is the subset of ref_lis
    """
    Determine the precursor mass. base_lis: theoretical intensity list, ref_lis: expanded experimental list 2 times window length,
    exp_lis: real peaks that are fragmented. rule1: final mass must include at least one non-mono iso-peak,
    rule2: only peaks within the mass window to be considered, rule3: assume the charge is correct,
    rule4: provided precursor mass is only used to derive the mass region, rule5: exp_lis only use the current MS1,
    rule6: peaks do not overlap, rule7: fuse_peak function makes sure only one fused peak in one possible region.
    return top 3 coefficient and calibrated mass
    :param base_lis:
    :param ref_lis:
    :param exp_lis:
    :param preCharge:
    :param ms1tol:
    """
    res = [(-1, 0), (-1, 0), (-1, 0)]  # large to small
    if not exp_lis:
        return res  # if nothing in the MS1 region, then nothing to be detect
    minShift = 0.99 / preCharge  # min mass shift from C, H, O , N, S
    maxShift = 1.0041 / preCharge  # max mass shift from C, H, O , N, S
    tolMz = 2 * ms1tol * exp_lis[0][0]  # m/z tolerance in dalton
    protonMass = 1.007276
    ave_unit = 1.0032  # average isotope mass in nature from human proteins
    lowerBound = exp_lis[0][0]
    upperBound = exp_lis[-1][0]
    left_ref = [peak for peak in ref_lis if peak[0] < exp_lis[0][0]]
    right_ref = [peak for peak in ref_lis if peak[0] >= exp_lis[0][0]]
    while exp_lis:
        curDeque = deque([exp_lis[0]])
        del exp_lis[0], right_ref[0]

        '''add left'''
        index = -1  # search from the end of left_ref
        while len(left_ref) >= -index and curDeque[0][0] - left_ref[index][0] <= maxShift + tolMz:
            if curDeque[0][0] - left_ref[index][0] >= minShift - tolMz:
                curDeque.appendleft(left_ref[index])
                del left_ref[index]
            else:
                index -= 1

        '''add right'''
        index = 0  # search from the begin of right_ref
        while len(right_ref) > index and right_ref[index][0] - curDeque[-1][0] <= maxShift + tolMz:
            if right_ref[index][0] - curDeque[-1][0] >= minShift - tolMz:
                curDeque.append(right_ref[index])
                if right_ref[index] in exp_lis:
                    exp_lis.remove(right_ref[index])
                del right_ref[index]
            else:
                index += 1

        intensity_list = [peak[1] for peak in curDeque]
        mz_list = [peak[0] for peak in curDeque]
        padding_intensity_list = [0] * len(base_lis) + intensity_list + [0] * (len(base_lis) - 1)
        start = bisect_left(mz_list, lowerBound)
        end = bisect_right(mz_list, upperBound)
        end_search = end + len(base_lis) - 2  # make sure non-mono inside
        for pos in range(start, end_search):
            curr_exp = padding_intensity_list[pos: pos + len(base_lis) + 1]
            curr_coef = np.corrcoef(curr_exp, [0] + base_lis)[0][1]

            if curr_coef > res[0][0]:
                weighted_left = max(pos - len(base_lis) + 1, 0)
                weighted_right = min(len(curDeque) - 1, pos)
                w_int = np.array(intensity_list[weighted_left: weighted_right + 1])
                w_mz = np.array(mz_list[weighted_left: weighted_right + 1]) - ave_unit / preCharge * \
                       np.arange(weighted_left + len(base_lis) - 1 - pos, weighted_right + len(base_lis) - pos)
                calibrated_mass = (np.average(w_mz, weights=w_int) - protonMass) * preCharge
                res[1] = res[0]
                res[2] = res[1]
                res[0] = (curr_coef, calibrated_mass)
            elif curr_coef > res[1][0]:
                weighted_left = max(pos - len(base_lis) + 1, 0)
                weighted_right = min(len(curDeque) - 1, pos)
                w_int = np.array(intensity_list[weighted_left: weighted_right + 1])
                w_mz = np.array(mz_list[weighted_left: weighted_right + 1]) - ave_unit / preCharge * \
                       np.arange(weighted_left + len(base_lis) - 1 - pos, weighted_right + len(base_lis) - pos)
                calibrated_mass = (np.average(w_mz, weights=w_int) - protonMass) * preCharge
                res[2] = res[1]
                res[1] = (curr_coef, calibrated_mass)
            elif curr_coef > res[2][0]:
                weighted_left = max(pos - len(base_lis) + 1, 0)
                weighted_right = min(len(curDeque) - 1, pos)
                w_int = np.array(intensity_list[weighted_left: weighted_right + 1])
                w_mz = np.array(mz_list[weighted_left: weighted_right + 1]) - ave_unit / preCharge * \
                       np.arange(weighted_left + len(base_lis) - 1 - pos, weighted_right + len(base_lis) - pos)
                calibrated_mass = (np.average(w_mz, weights=w_int) - protonMass) * preCharge
                res[2] = (curr_coef, calibrated_mass)

    return res


class Ms2PrecursorInfo:
    __slots__ = 'preCharge', 'preMz', 'currentNum', 'preNum', 'preIntensity'

    def __init__(self, pre_charge, pre_mz, curr_num, pre_num):
        self.preCharge = pre_charge
        self.preMz = pre_mz
        self.currentNum = curr_num
        self.preNum = pre_num


def fuse_peaks(peaks, ms1tol, preCharge):
    """
    Merge peaks if they are within the isotope mass tolerance
    :param peaks:
    :param ms1tol:
    :param preCharge:
    """
    diff = (1.0041 - 0.99) / preCharge  # max isotopic mass difference for C, H, O, N, S
    peaks = sorted(peaks)
    res = peaks[0:1]
    start = 1
    while start < len(peaks):
        if np.isclose(peaks[start][0], res[-1][0], rtol=2 * ms1tol, atol=diff):
            tempX = (peaks[start][0] * peaks[start][1] + res[-1][0] * res[-1][1]) / \
                    (peaks[start][1] + res[-1][1])
            tempY = peaks[start][1] + res[-1][1]
            res[-1] = (tempX, tempY)
        else:
            res.append(peaks[start])
        start += 1
    return res


def pre_mass_refine_average_union(tuple_bags):
    """
    Return the dictionary of refined precursor mass{'scan_num1':mass1, 'scan_num2': mass2, ...}
    :param tuple_bags:
    """
    path, ms1tol, coe_threshold = tuple_bags
    iso_win = 2.0  # Thomson unit, isolation window for precursor ion extraction

    with mzxml.read(path) as spectra:

        ms1dict = dict()  # create ms1 dictionary {'num':[(mz1,int1),(ms2,int2),(mz3,int3)]}
        ms1Linked_list = []  # create MS1 linked list
        ms2Info = []  # store each MS2 information
        for spectrum in spectra:
            if spectrum['msLevel'] == 1:
                ms1dict[int(spectrum['num'])] = [spectrum['m/z array'], spectrum['intensity array']]
                ms1Linked_list.append(int(spectrum['num']))
            elif spectrum['msLevel'] == 2:

                '''store the ms2 information'''
                ms2Info.append(Ms2PrecursorInfo(spectrum['precursorMz'][0]['precursorCharge'],
                                                spectrum['precursorMz'][0]['precursorMz'],
                                                spectrum['num'], spectrum['precursorMz'][0]['precursorScanNum']))

    '''start to calculate the correct precursor mass, current ms -> previous ms -> latter ms -> unchanged
    the minimum mass shift is 0.99 and maximum is 1.0041 (from Kojak) and the average mass is 1.0032 (from pParse).
    the isotopic cluster should be continuously extracted, but the mono mass can be missed'''
    res = dict()  # result dictionary
    lastCache = {'m/z': 0.0, 'charge': 0, 'preMass': 0.0}  # store the last result for later use, designed for HCD ETD situation
    for ms2 in ms2Info:
        if np.isclose(ms2.preMz, lastCache['m/z']) and ms2.preCharge == lastCache['charge']:
            res[ms2.currentNum] = lastCache['preMass']
            continue
        theoIsotopeCluster, maxRatio = binomial(ms2.preMz * ms2.preCharge)
        '''find the ms1 and the labeled peak, limit the peak range in (m/z-redundancy, m/z + redundancy + maximumNum)'''
        preScanIdxLeft = bisect_left(ms1Linked_list, int(ms2.preNum))
        preScanIdxRight = preScanIdxLeft + 1
        retentionScan = 30
        maxNum = 4
        fusePeaks = []
        '''add ms1 from left'''
        while preScanIdxLeft >= 0 and int(ms2.preNum) - ms1Linked_list[preScanIdxLeft] < retentionScan and maxNum > 0:
            refWindow = pick_peak(ms1dict[ms1Linked_list[preScanIdxLeft]], ms2.preMz, 2 * iso_win)
            fusePeaks += refWindow
            preScanIdxLeft -= 1
            maxNum -= 1
        maxNum = 3
        '''add ms1 from right'''
        while preScanIdxRight < len(ms1Linked_list) and \
                ms1Linked_list[preScanIdxRight] - int(ms2.preNum) < retentionScan and maxNum > 0:
            refWindow = pick_peak(ms1dict[ms1Linked_list[preScanIdxRight]], ms2.preMz, 2 * iso_win)
            fusePeaks += refWindow
            preScanIdxRight += 1
            maxNum -= 1

        refWindow = fuse_peaks(fusePeaks, ms1tol, ms2.preCharge)
        expWindow = [peak for peak in refWindow if
                     ms2.preMz + iso_win / 2 >= peak[0] >= ms2.preMz - iso_win / 2]

        finalRes = precursor_determination(theoIsotopeCluster, refWindow, expWindow, ms2.preCharge, ms1tol)
        finalRes = top_res_merge([], finalRes, ms1tol)
        finalRes = [finalRes[0]] + [fe for fe in finalRes[1:] if fe[0] > coe_threshold]

        if finalRes[0][0] < coe_threshold:
            finalRes = [(1, (ms2.preMz - 1.00728) * ms2.preCharge)]

        finalRes = [ii[1] for ii in finalRes]
        res[ms2.currentNum] = finalRes
        lastCache['preMass'] = res[ms2.currentNum]
        lastCache['charge'] = ms2.preCharge
        lastCache['m/z'] = ms2.preMz

    with open(r'pickled_pair_wise/{}'.format(path.split('/')[-1].split('.')[0]), 'rb') as file:
        pair_wise_data = pickle.load(file)
    new_pair_wise_data = []
    for ele in pair_wise_data:
        if str(ele['CID_scan']) in res.keys():
            if len(res[str(ele['CID_scan'])]) > 1:
                for idx in range(1, len(res[str(ele['CID_scan'])]) + 1):
                    new_ele = {'mass': res[str(ele['CID_scan'])][idx - 1], 'CID_scan': ele['CID_scan'] + idx * 0.1,
                               'charge': ele['charge'], 'CID_peaks': ele['CID_peaks'], 'ETD_peaks': ele['ETD_peaks']}
                    new_pair_wise_data.append(new_ele)
            else:
                ele['mass'] = res[str(ele['CID_scan'])][0]
                new_pair_wise_data.append(ele)

    with open(r'pickled_pair_wise/{}'.format(path.split('/')[-1].split('.')[0]), 'wb') as file:
        pickle.dump(new_pair_wise_data, file)
    return res


if __name__ == "__main__":
    """Highly recommend to run this module right after the spectra separation to refine the precursor mass.
     although it is not compulsory to do so."""

    with open('ECLPF_conf', 'r') as conf_file:
        conf = json.load(conf_file)
    TOLERANCE1 = conf['ms1_tol']
    for ele in conf['data_path']:
        if ' ' in ele:
            print('space contains')
            quit()
    threshold = 0.9  # pearson correlation threshold to trust the refined result
    os.chdir('data/')
    args = []  # for multi-processing
    for i in conf['data_path']:
        args.append((i, TOLERANCE1, threshold))

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        result = executor.map(pre_mass_refine_average_union, args)
