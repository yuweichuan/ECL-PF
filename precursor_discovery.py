import pickle
import numpy as np
import time
import os
import concurrent.futures
import shutil


def find_precursor_ab_mass(sorted_list, pre_mass, xl, ml, ms, tol, signum=1):
    # sorted_list = [(mass1,intensity1),(ma2,int2),()], tol unit is ms2 ppm 2e-5
    """return list of tuples of potential alpha and beta precursor masses,
    output = [[beta,alpha,#markers,intensity],[],[]],  sort from higher possible ones to lower ones"""
    if len(sorted_list) < 2:
        return [[0, 0, -1, 0]]  # return 0 result
    delta = ml - ms
    res = []
    indx1 = 0
    indx2 = 1
    marker = True
    backup_sorted_list = list(sorted_list)
    while marker:
        if sorted_list[indx2][0] - sorted_list[indx1][0] > \
                delta + tol * (sorted_list[indx2][0] + sorted_list[indx1][0]):
            indx1 += 1
        elif sorted_list[indx2][0] - sorted_list[indx1][0] < \
                delta - tol * (sorted_list[indx2][0] + sorted_list[indx1][0]):
            indx2 += 1
        else:
            weighted_mass = ((sorted_list[indx1][0] - ms) * sorted_list[indx1][1] +
                             (sorted_list[indx2][0] - ml) * sorted_list[indx2][1]) /\
                            (sorted_list[indx1][1] + sorted_list[indx2][1])
            res.append((round(weighted_mass, 3), sorted_list[indx1][1] + sorted_list[indx2][1]))
            # free mass of peptide and summation of their marker intensities
            if sorted_list[indx1] in backup_sorted_list:
                backup_sorted_list.remove(sorted_list[indx1])
            backup_sorted_list.remove(sorted_list[indx2])
            # remove these markers for the following usage
            indx1 += 1
            indx2 += 1
        if indx2 >= len(sorted_list):
            marker = False
    validate_res = []
    if len(res) >= 1:
        first = 0
        last = len(res) - 1
        while first <= last:
            if res[first][0] + res[last][0] < pre_mass - xl - tol * (res[first][0] + res[last][0]):
                first += 1
            elif res[first][0] + res[last][0] > pre_mass - xl + tol * (res[first][0] + res[last][0]):
                last -= 1
            else:
                validate_res.append([res[first][0], res[last][0], 4, res[first][1] + res[last][1]])
                # type_four adding
                last -= 1
                first += 1
        if validate_res:  # remove marker peaks if 4 marker found
            for val_res_ele in validate_res:
                rem1 = [i for i in res if i[0] == val_res_ele[0]]
                rem2 = [i for i in res if i[0] == val_res_ele[1]]
                if rem1 != rem2:
                    res.remove(rem1[0])
                    res.remove(rem2[0])
                else:
                    res.remove(rem1[0])

        for res_ele in res:
            sorted_list_mz = [i[0] for i in sorted_list]
            req1 = np.isclose(pre_mass - xl - res_ele[0] + ml, sorted_list_mz, rtol=tol)
            req2 = np.isclose(pre_mass - xl - res_ele[0] + ms, sorted_list_mz, rtol=tol)
            req = req1 + req2
            # check if there exist the third solo peak
            if any(req):
                indx_3_peak = [peak3 for peak3 in range(len(req)) if req[peak3]]
                for index3 in indx_3_peak:
                    if sorted_list[index3] in backup_sorted_list:
                        backup_sorted_list.remove(sorted_list[index3])  # remove the third solo peak
                        validate_res.append([res_ele[0], round(pre_mass - xl - res_ele[0], 3),
                                             3, res_ele[1] + sorted_list[index3][1]])
                    else:
                        validate_res.append([res_ele[0], round(pre_mass - xl - res_ele[0], 3),
                                             3.5, res_ele[1] + sorted_list[index3][1]])  # define the extreme 3.5 type
                # type_three adding
            else:
                validate_res.append([res_ele[0], round(pre_mass - xl - res_ele[0], 3), 2, res_ele[1]])
                # type_two adding
    backup_type_one = list(backup_sorted_list)  # for the use of type one
    """start to find 1.25 peaks in one from longer arm and the other from shorter arm"""
    indx_start = 0
    indx_end = len(backup_sorted_list) - 1
    while indx_start < indx_end:
        if backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0] < \
                pre_mass - tol * (backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0]):
            indx_start += 1
        elif backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0] > \
                pre_mass + tol * (backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0]):
            indx_end -= 1
        else:
            validate_res.append([round(backup_sorted_list[indx_start][0] - ms, 3),
                                 round(backup_sorted_list[indx_end][0] - ml, 3),
                                 1.25, backup_sorted_list[indx_start][1] + backup_sorted_list[indx_end][1]])
            validate_res.append([round(backup_sorted_list[indx_start][0] - ml, 3),
                                 round(backup_sorted_list[indx_end][0] - ms, 3),
                                 1.25, backup_sorted_list[indx_start][1] + backup_sorted_list[indx_end][1]])
            # adding the first situation type_1.5
            backup_type_one.remove(backup_sorted_list[indx_start])
            backup_type_one.remove(backup_sorted_list[indx_end])
            indx_start += 1
            indx_end -= 1
    """start to find 1.5 peaks in both from shorter arm"""
    indx_start = 0
    indx_end = len(backup_sorted_list) - 1
    while indx_start <= indx_end:
        if backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0] + delta < \
                pre_mass - tol * (backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0]):
            indx_start += 1
        elif backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0] + delta > \
                pre_mass + tol * (backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0]):
            indx_end -= 1
        else:
            validate_res.append([round(backup_sorted_list[indx_start][0] - ms, 3),
                                 round(backup_sorted_list[indx_end][0] - ms, 3),
                                 1.5, backup_sorted_list[indx_start][1] + backup_sorted_list[indx_end][1]])

            # adding the third situation type_1.5
            backup_type_one = [i for i in backup_type_one if i not in
                               [backup_sorted_list[indx_start], backup_sorted_list[indx_end]]]
            indx_start += 1
            indx_end -= 1
    """start to find 1.5 peaks in both from longer arm"""
    indx_start = 0
    indx_end = len(backup_sorted_list) - 1
    while indx_start <= indx_end:
        if backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0] - delta < \
                pre_mass - tol * (backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0]):
            indx_start += 1
        elif backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0] - delta > \
                pre_mass + tol * (backup_sorted_list[indx_start][0] + backup_sorted_list[indx_end][0]):
            indx_end -= 1
        else:
            validate_res.append([round(backup_sorted_list[indx_start][0] - ml, 3),
                                 round(backup_sorted_list[indx_end][0] - ml, 3),
                                 1.5, backup_sorted_list[indx_start][1] + backup_sorted_list[indx_end][1]])

            # adding the second situation type_1.5
            backup_type_one = [i for i in backup_type_one if i not in
                               [backup_sorted_list[indx_start], backup_sorted_list[indx_end]]]
            indx_start += 1
            indx_end -= 1
    """start to find type one peak"""
    for peak1 in backup_type_one:
        mass_sit1 = round(peak1[0] - ms, 3)  # situation 1, peak from shorter arm
        mass_sit2 = round(peak1[0] - ml, 3)  # situation 2, peak from longer arm
        validate_res.append([mass_sit1, round(pre_mass - xl - mass_sit1, 3), 1, peak1[1]])
        validate_res.append([mass_sit2, round(pre_mass - xl - mass_sit2, 3), 1, peak1[1]])
    """sort from beta to alpha"""
    for val_res_ord in range(len(validate_res)):
        if validate_res[val_res_ord][0] > validate_res[val_res_ord][1]:
            validate_res[val_res_ord] = [validate_res[val_res_ord][1], validate_res[val_res_ord][0],
                                         validate_res[val_res_ord][2], validate_res[val_res_ord][3]]
    """mass could not be negative ones """
    validate_res = [i for i in validate_res if i[0] > 0]
    """sort by number of markers, then followed by intensity"""
    validate_res = sorted(validate_res, key=lambda x: (-x[2], -x[3]))
    """only return high possible ones if exist markers >= 3, else return everything"""
    # if validate_res and validate_res[0][2] >= 3:
    #     validate_res = [val_res_ele for val_res_ele in validate_res if val_res_ele[2] >= 3]
    validate_res = [val_res_ele for val_res_ele in validate_res if val_res_ele[2] >= min(signum, 4)]
    return validate_res  # [[beta,alpha,#markers,intensity],[],[]],  sort from higher possible ones to lower ones


def extract_peptide(chain_mass, tol_da, path='database_file/'):  # mass precision in 4 decimal
    """given a concrete mass and tolerance,
    return possible peptides in the database [(mass,'PEPTIDE',description),( , , ),( , , )...]"""
    dir_mass = os.listdir(path)
    mass_idx = [(chain_mass - int(dir_mass[i]), i) for i in range(len(dir_mass)) if chain_mass - int(dir_mass[i]) >= 0]
    with open(path + dir_mass[min(mass_idx)[1]], 'rb') as file:
        mass_list = pickle.load(file)
    pep_keys = set(np.around(np.arange(chain_mass - tol_da, chain_mass + tol_da, 0.001), 3))
    pep_db = {key: mass_list[1][key] for key in mass_list[0] & pep_keys}
    res = []
    for key, values in pep_db.items():
        for value in values:
            res.append((key, value[0], value[1]))
    return res


def db_to_spectra_extraction(correspond_matrix, tol, path='database_file/'):
    """from database to extract chain precursor candidates
    returns list [[[beta,alpha,#,int],[beta,alpha,#,int],[],..],
                        [[beta,alpha,#,int],[beta,alpha,#,int],[],..]
                        ...]
                        where beta = [(mass,concatenate_sequence,concatenate_des),(mass,con_seq,con_des),(),...]
                        alpha = [(mass,concatenate_sequence,concatenate_des),(mass,con_seq,con_des),(),...]"""

    dir_mass = os.listdir(path)
    num_dir_mass = sorted([int(dir_mass_ele) for dir_mass_ele in dir_mass])
    mass_interval = num_dir_mass[1] - num_dir_mass[0]
    """create target list"""
    res = correspond_matrix
    """replace the mass element with candidates in res"""
    for ele in dir_mass:
        with open(path + ele, 'rb') as file:
            sub_mass = pickle.load(file)  # [set(), dict()]
        for key in res:
            for key_ele in key:
                if not isinstance(key_ele[0], list) and mass_interval > key_ele[0] - float(ele) >= 0:
                    beta_list = []
                    candidate_set = set(np.round(np.arange(key_ele[0] * (1 - tol), key_ele[0] * (1 + tol), 0.001), 3))
                    candidate_key = candidate_set & sub_mass[0]
                    # could be empty key
                    for candidate_key_ele in candidate_key:
                        beta_list.append(
                            (candidate_key_ele, sub_mass[1][candidate_key_ele][0], sub_mass[1][candidate_key_ele][1]))
                    key_ele[0] = beta_list

                if not isinstance(key_ele[1], list) and mass_interval > key_ele[1] - float(ele) >= 0:
                    alpha_list = []
                    candidate_set = set(np.round(np.arange(key_ele[1] * (1 - tol), key_ele[1] * (1 + tol), 0.001), 3))
                    candidate_key = candidate_set & sub_mass[0]
                    # could be empty key
                    for candidate_key_ele in candidate_key:
                        alpha_list.append(
                            (candidate_key_ele, sub_mass[1][candidate_key_ele][0], sub_mass[1][candidate_key_ele][1]))
                    key_ele[1] = alpha_list
    """replace the out of range mass with empty list"""
    for key in res:
        for key_ele in key:
            if not isinstance(key_ele[0], list):
                key_ele[0] = []
            if not isinstance(key_ele[1], list):
                key_ele[1] = []

    return res


if __name__ == '__main__':
    with open(r"F:\OneDrive - HKUST Connect\Research\ECL_X\ECLX_src\data\pickled_pair_wise\MS200295-CBDPS-HCD",'rb') as f:
        spectra = pickle.load(f)
    print(f'spcetra number is {len(spectra)} \n the first element is:')
    print(spectra[0:1])
    xl_mass = 509.0963
    large_mass = 455.0868
    small_mass = 54.0106
    TOL = 2e-5
    pre_chian_tol = 3e-5








