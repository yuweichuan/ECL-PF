import numpy as np
from pyteomics import mass
import concurrent.futures
import os
import json
import pickle
import time
from precursor_discovery import find_precursor_ab_mass, db_to_spectra_extraction
from match import spec_pre_process, cid_x_corr_score, cid_etd_x_corr_score, \
    cid_p_score, cid_etd_p_score, cid_merox_score, cid_etd_merox_score
from fragment import fast_b_y_fragment, fast_b_y_c_z_fragment
import csv
from statistics import *
from fdr import fdr
from protein_score import protein_score_filtering, export_all
import logging


class TopCandidate:
    __slots__ = 'alphaSet', 'betaSet', 'marker'

    def __init__(self, alpha1, beta1, marker):
        self.alphaSet = alpha1
        self.betaSet = beta1
        self.marker = marker


def splitPossibleResult(topAlpha, topBeta, marker, scan, pre_mass):
    """split the possible results into individual ones, a2, b2 is the average score here,
    return list of dicts: [{'CID_scan': 123123, 'mass': 799, alpha1 : [-1, 'PEPTIDE', -1, ['0'], 799],
    beta1=[...], alpha2=a2 score, beta2= b2 score},{}]"""
    res = []
    flag = []
    a2 = []
    b2 = []
    for _alpha in topAlpha:
        for _beta in topBeta:
            a2.append(_alpha[0])
            b2.append(_beta[0])
            if (_alpha[1], _beta[1]) not in flag:
                flag.append((_alpha[1], _beta[1]))
                res.append({'CID_scan': scan, 'mass': pre_mass, 'marker': marker, 'alpha1': list(_alpha),
                            'beta1': list(_beta)})
            else:
                index_ab = flag.index((_alpha[1], _beta[1]))
                res[index_ab]['alpha1'][3].extend(_alpha[3])
                res[index_ab]['beta1'][3].extend(_beta[3])
                res[index_ab]['alpha1'][3] = list(set(res[index_ab]['alpha1'][3]))
                res[index_ab]['beta1'][3] = list(set(res[index_ab]['beta1'][3]))
    a2 = sum(a2) / len(a2)
    b2 = sum(b2) / len(b2)
    for _res in res:
        _res['alpha2'] = a2
        _res['beta2'] = b2

    return res


def cid_etd_overlap(sequence, alphaPos, alphaMass, preMass, massDict, tol, cidSpec, etdSpec):  # ms2 tol ppm
    """return lists (cid, etd) of tuple of the matched ions [(mass,intensity),(,)...]
    and the reduced spectrum (cid, etd)"""
    _newCidSpec = list(cidSpec)
    _newEtdSpec = list(etdSpec)
    _cidMatchedList = []
    _etdMatchedList = []
    cid_theo_list, etd_theo_list = \
        fast_b_y_c_z_fragment(sequence, alphaPos, alphaMass, preMass, massDict)

    index_theo = 0
    index_exp = 0
    while index_exp < len(cidSpec) and index_theo < len(cid_theo_list):
        if cidSpec[index_exp][0] - cid_theo_list[index_theo] < \
                - tol * max(cidSpec[index_exp][0], cid_theo_list[index_theo]):
            index_exp += 1
        elif cidSpec[index_exp][0] - cid_theo_list[index_theo] > \
                tol * max(cidSpec[index_exp][0], cid_theo_list[index_theo]):
            index_theo += 1
        else:
            _cidMatchedList.append(cidSpec[index_exp])
            _newCidSpec[index_exp] = (_newCidSpec[index_exp][0], _newCidSpec[index_exp][1] / 2)  # divide intensity by 2
            index_exp += 1
            index_theo += 1

    index_theo = 0
    index_exp = 0
    while index_exp < len(etdSpec) and index_theo < len(etd_theo_list):
        if etdSpec[index_exp][0] - etd_theo_list[index_theo] < \
                - tol * max(etdSpec[index_exp][0], etd_theo_list[index_theo]):
            index_exp += 1
        elif etdSpec[index_exp][0] - etd_theo_list[index_theo] > \
                tol * max(etdSpec[index_exp][0], etd_theo_list[index_theo]):
            index_theo += 1
        else:
            _etdMatchedList.append(etdSpec[index_exp])
            _newEtdSpec[index_exp] = (_newEtdSpec[index_exp][0], _newEtdSpec[index_exp][1] / 2)  # divide intensity by 2
            index_exp += 1
            index_theo += 1

    return _cidMatchedList, _etdMatchedList, _newCidSpec, _newEtdSpec


def cid_overlap(sequence, alphaPos, alphaMass, preMass, massDict, tol, cidSpec):  # ms2 tol ppm
    """return lists (cid) of tuple of the matched ions [(mass,intensity),(,)...]
    and the reduced spectrum (cid)"""
    _newCidSpec = list(cidSpec)
    _cidMatchedList = []
    cid_theo_list = fast_b_y_fragment(sequence, alphaPos, alphaMass, preMass, massDict)

    index_theo = 0
    index_exp = 0
    while index_exp < len(cidSpec) and index_theo < len(cid_theo_list):
        if cidSpec[index_exp][0] - cid_theo_list[index_theo] < \
                - tol * max(cidSpec[index_exp][0], cid_theo_list[index_theo]):
            index_exp += 1
        elif cidSpec[index_exp][0] - cid_theo_list[index_theo] > \
                tol * max(cidSpec[index_exp][0], cid_theo_list[index_theo]):
            index_theo += 1
        else:
            _cidMatchedList.append(cidSpec[index_exp])
            _newCidSpec[index_exp] = (_newCidSpec[index_exp][0], _newCidSpec[index_exp][1] / 2)  # divide intensity by 2
            index_exp += 1
            index_theo += 1

    return _cidMatchedList, _newCidSpec


def insertion(candidate, topCandidates, topT):
    """candidate format = (score, sequence, xl position, description, candidate mass),
    topCandidates ordered by the score from smaller to larger, insert candidate to the proper position,
    return the topCandidates with same length if possible"""
    if candidate[0] < topCandidates[0][0]:  # we need to consider tie to the first place situation
        return topCandidates
    elif candidate[0] >= topCandidates[-1][0]:
        topCandidates.append(candidate)
        if topCandidates[-topT][0] != topCandidates[-1 - topT][0]:
            return topCandidates[-topT:]
        else:
            return topCandidates
    else:
        lo = 0
        hi = len(topCandidates)
        while lo < hi:
            mid = (lo+hi)//2
            if topCandidates[mid][0] < candidate[0]:
                lo = mid+1
            else:
                hi = mid
        topCandidates.insert(lo, candidate)
        if topCandidates[-topT][0] != topCandidates[-1 - topT][0]:
            return topCandidates[-topT:]
        else:
            return topCandidates


def run(tuple_bag):
    sub_spectra, tol1, tol2, xl, ml, ms, link_site, mass_dict = tuple_bag
    '''find alpha beta chain mass of sub spectra
    returns [{'CID_scan':_, 'mass':_, 'alpha1':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'beta1':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'alpha2':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'beta2':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'marker':_,},{},{}]'''
    t_alpha = 30  # top number of alpha candidates
    t_beta = 30  # top number of beta candidates
    truncation = 0.6  # top numbers' scores should be in a range, cannot be too different
    res = []
    correspond_matrix = []

    for ele in sub_spectra:
        if 'rep_CID_peaks' in ele.keys():
            sorted_list = [peak for peak in ele['rep_CID_peaks'] if peak[0] < ele['mass']]
            correspond_matrix.append(find_precursor_ab_mass(sorted_list, ele['mass'], xl, ml, ms, tol2,
                                                            signum=ele['repNum']))
        else:
            sorted_list = [peak for peak in ele['CID_peaks'] if peak[0] < ele['mass']]
            correspond_matrix.append(find_precursor_ab_mass(sorted_list, ele['mass'], xl, ml, ms, tol2, signum=1))

    candidate_list = db_to_spectra_extraction(correspond_matrix, tol1 + tol2, path='database_file/')
    for idx in range(len(sub_spectra)):  # calculate every spectrum
        ele = sub_spectra[idx]
        preprocessIn = [peak for peak in ele['CID_peaks'] if peak[0] < ele['mass']]
        if not preprocessIn:
            continue
        process_CID_spec = spec_pre_process(preprocessIn)
        process_ETD_spec = []  # to be calculated
        etdConsider = False
        if ele['ETD_peaks']:
            process_ETD_spec = spec_pre_process(ele['ETD_peaks'])
            etdConsider = True

        topAlpha = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_alpha
        topBeta = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_beta
        marker = 0
        # hit(s) with score,sequence,position,[description],mass

        if etdConsider:
            for CID_sub_candidate in candidate_list[idx]:  # [beta,alpha,#,int], 3rd layer loop

                if CID_sub_candidate[0] and CID_sub_candidate[1]:
                    specMassRegionTopAlpha = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_alpha
                    specMassRegionTopBeta = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_beta

                    for con_alpha in CID_sub_candidate[1]:  # concatenated alpha sequence and description
                        alpha_seq_uncompress = con_alpha[1].split('|')
                        alpha_des_uncompress = con_alpha[2].split('|')
                        alpha_num = len(alpha_seq_uncompress)
                        for idx_con in range(alpha_num):  # 1st layer loop
                            alpha = (con_alpha[0], alpha_seq_uncompress[idx_con], alpha_des_uncompress[idx_con])
                            # mass, sequence, description
                            '''xcorr score'''
                            # alpha_candidate = cid_etd_x_corr_score(process_CID_spec, process_ETD_spec, alpha[1],
                            #                                        ele['mass'], alpha[0], tol2, mass_dict, link_site, proN=True)
                            '''p score in xlinkx'''
                            alpha_candidate = cid_etd_p_score(process_CID_spec, process_ETD_spec, alpha[1], ele['mass'],
                                                              alpha[0], tol2, mass_dict, link_site, xl, proN=True)
                            '''MeroX score'''
                            # alpha_candidate = cid_etd_merox_score(ele['CID_peaks'], ele['ETD_peaks'], alpha[1],
                            #                                       ele['mass'], alpha[0], tol2, mass_dict, link_site,
                            #                                       ml, ms, ele['charge'], proN=True)

                            # (score, sequence, xl position)
                            if alpha_candidate:
                                currAlpha = (alpha_candidate[0], alpha_candidate[1], alpha_candidate[2], alpha[2], alpha[0])
                                specMassRegionTopAlpha = insertion(currAlpha, specMassRegionTopAlpha, t_alpha)

                    for con_beta in CID_sub_candidate[0]:  # concatenate beta sequence and description, 2nd layer loop
                        beta_seq_uncompress = con_beta[1].split('|')
                        beta_des_uncompress = con_beta[2].split('|')
                        beta_num = len(beta_seq_uncompress)
                        for idx_con in range(beta_num):  # 1st layer loop
                            beta = (con_beta[0], beta_seq_uncompress[idx_con], beta_des_uncompress[idx_con])
                            # mass, sequence, description
                            '''xcorr score'''
                            # beta_candidate = cid_etd_x_corr_score(process_CID_spec, process_ETD_spec, beta[1],
                            #                                       ele['mass'], beta[0], tol2, mass_dict, link_site, proN=True)
                            '''p score in xlinkx'''
                            beta_candidate = cid_etd_p_score(process_CID_spec, process_ETD_spec, beta[1], ele['mass'],
                                                             beta[0], tol2, mass_dict, link_site, xl, proN=True)
                            '''MeroX score'''
                            # beta_candidate = cid_etd_merox_score(ele['CID_peaks'], ele['ETD_peaks'], beta[1],
                            #                                      ele['mass'], beta[0], tol2, mass_dict, link_site,
                            #                                      ml, ms, ele['charge'], proN=True)

                            # (score, sequence, xl position)
                            if beta_candidate:
                                currBeta = (beta_candidate[0], beta_candidate[1], beta_candidate[2], beta[2], beta[0])
                                specMassRegionTopBeta = insertion(currBeta, specMassRegionTopBeta, t_beta)

                    '''compare the score with the highest'''
                    if specMassRegionTopBeta[-1][0] + specMassRegionTopAlpha[-1][0] > topAlpha[-1][0] + topBeta[-1][0]:
                        topAlpha = specMassRegionTopAlpha
                        topBeta = specMassRegionTopBeta
                        marker = CID_sub_candidate[2]

        else:
            for CID_sub_candidate in candidate_list[idx]:  # [beta,alpha,#,int], 3rd layer loop

                if CID_sub_candidate[0] and CID_sub_candidate[1]:
                    specMassRegionTopAlpha = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_alpha
                    specMassRegionTopBeta = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_beta

                    for con_alpha in CID_sub_candidate[1]:  # concatenated alpha sequence and description
                        alpha_seq_uncompress = con_alpha[1].split('|')
                        alpha_des_uncompress = con_alpha[2].split('|')
                        alpha_num = len(alpha_seq_uncompress)
                        for idx_con in range(alpha_num):  # 1st layer loop
                            alpha = (con_alpha[0], alpha_seq_uncompress[idx_con], alpha_des_uncompress[idx_con])
                            # mass, sequence, description
                            '''xcorr score'''
                            # alpha_candidate = cid_x_corr_score(process_CID_spec, alpha[1],
                            #                                    ele['mass'], alpha[0], tol2, mass_dict, link_site, proN=True)
                            '''p score in xlinkx'''
                            alpha_candidate = cid_p_score(process_CID_spec, alpha[1], ele['mass'], alpha[0], tol2,
                                                          mass_dict, link_site, xl, proN=True)
                            '''MeroX score'''
                            # alpha_candidate = cid_merox_score(ele['CID_peaks'], alpha[1], ele['mass'], alpha[0], tol2,
                            #                                   mass_dict, link_site, ml, ms, ele['charge'], proN=True)


                            # (score, sequence, xl position)
                            if alpha_candidate:
                                currAlpha = (alpha_candidate[0], alpha_candidate[1], alpha_candidate[2], alpha[2], alpha[0])
                                specMassRegionTopAlpha = insertion(currAlpha, specMassRegionTopAlpha, t_alpha)

                    for con_beta in CID_sub_candidate[0]:  # concatenate beta sequence and description, 2nd layer loop
                        beta_seq_uncompress = con_beta[1].split('|')
                        beta_des_uncompress = con_beta[2].split('|')
                        beta_num = len(beta_seq_uncompress)
                        for idx_con in range(beta_num):  # 1st layer loop
                            beta = (con_beta[0], beta_seq_uncompress[idx_con], beta_des_uncompress[idx_con])
                            # mass, sequence, description
                            '''xcorr score'''
                            # beta_candidate = cid_x_corr_score(process_CID_spec, beta[1],
                            #                                   ele['mass'], beta[0], tol2, mass_dict, link_site, proN=True)
                            '''p score in xlinkx'''
                            beta_candidate = cid_p_score(process_CID_spec, beta[1], ele['mass'], beta[0], tol2,
                                                         mass_dict, link_site, xl, proN=True)
                            '''MeroX score'''
                            # beta_candidate = cid_merox_score(ele['CID_peaks'], beta[1], ele['mass'], beta[0], tol2,
                            #                                  mass_dict, link_site, ml, ms, ele['charge'], proN=True)


                            # (score, sequence, xl position)
                            if beta_candidate:
                                currBeta = (beta_candidate[0], beta_candidate[1], beta_candidate[2], beta[2], beta[0])
                                specMassRegionTopBeta = insertion(currBeta, specMassRegionTopBeta, t_beta)

                    '''compare the score with the highest'''
                    if specMassRegionTopBeta[-1][0] + specMassRegionTopAlpha[-1][0] > topAlpha[-1][0] + topBeta[-1][0]:
                        topAlpha = specMassRegionTopAlpha
                        topBeta = specMassRegionTopBeta
                        marker = CID_sub_candidate[2]
        # splitRes = splitPossibleResult(topCandidates, alpha2, beta2, ele['CID_scan'], ele['mass'])
        # res.extend(splitRes)
        '''truncate the top candidates'''
        topAlpha = [_alpha for _alpha in topAlpha if np.isclose(_alpha[0], topAlpha[-1][0], atol=0, rtol=truncation)]
        topBeta = [_beta for _beta in topBeta if np.isclose(_beta[0], topBeta[-1][0], atol=0, rtol=truncation)]
        '''add to the result'''
        res.append((topAlpha, topBeta, marker, ele['CID_scan'], ele['mass']))

    return res


def all_abs_path_of_components(path=r'data\pickled_pair_wise'):
    res = []
    for maindir, subdir, file_name_list in os.walk(path):
        for filename in file_name_list:
            apath = os.path.join(maindir, filename)
            res.append(apath)

    return res


if __name__ == '__main__':
    '''input of the software'''
    with open('ECLPF_conf', 'r') as conf_file:
        conf = json.load(conf_file)

    Link_site = conf['link_site']  # support multiple sites such as ['K', 'R']
    Fix_mod = conf['fix_mod']
    Var_mod = conf['var_mod']
    MS1_TOLERANCE = conf['ms1_tol']
    MS2_TOLERANCE = conf['ms2_tol']
    xl_mass = conf['xl_mass']
    m_short = conf['m_short']
    m_long = conf['m_long']

    '''add modification masses into the dict'''
    new_mass = {}
    nTerm = {}
    cTerm = {}

    for i, j in Var_mod.items():
        new_mass[i] = j[0]
        if ('Peptide-nterm' in j[1]) or ('Protein-nterm' in j[1]):
            nTerm[i] = j[0]
        if ('Peptide-cterm' in j[1]) or ('Protein-cterm' in j[1]):
            cTerm[i] = j[0]

    for i, j in Fix_mod.items():
        new_mass[i] = j[0]
        for ele in nTerm:
            new_mass[ele + ')(' + i] = j[0] + nTerm[ele]
        for ele in cTerm:
            new_mass[i + ')(' + ele] = j[0] + nTerm[ele]

    for i, j in Var_mod.items():
        for ele in nTerm:
            new_mass[ele + ')(' + i] = j[0] + nTerm[ele]
        for ele in cTerm:
            new_mass[i + ')(' + ele] = j[0] + nTerm[ele]

    AA_mass = dict(mass.std_aa_mass)
    AA_mass.update(new_mass)

    '''import spectra file, after run spectra_separation'''

    spectra_list = all_abs_path_of_components()
    print(spectra_list)
    for spectra_path in spectra_list:
        with open(spectra_path, 'rb') as f:
            spectra = pickle.load(f)
        spectra_div = 500

        section = [i for i in range(len(spectra)) if i % spectra_div == 0]
        section_length = len(section)
        args = []
        for i in range(section_length):
            if i != len(section) - 1:
                args.append((spectra[section[i]: section[i + 1]], MS1_TOLERANCE, MS2_TOLERANCE, xl_mass,
                             m_long, m_short, Link_site, AA_mass))
            else:
                args.append((spectra[section[i]:], MS1_TOLERANCE, MS2_TOLERANCE, xl_mass,
                             m_long, m_short, Link_site, AA_mass))
        print('start multiprocessing')
        t1 = time.perf_counter()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            result = executor.map(run, args)

        t2 = time.perf_counter()
        time_diff = round(t2 - t1, 3)
        print(f"It took {time_diff} Secs to execute this method")
        t1 = time.perf_counter()
        with open(r'protein_name', 'rb') as f:
            protein_name = pickle.load(f)
        numberK = {key: protein_name[key][1] for key in protein_name.keys()}

        psm = protein_score_filtering(result, MS1_TOLERANCE, xl_mass, numberK)
        # psm = export_all(result)
        t2 = time.perf_counter()
        print('{} seconds to filter the result by protein score, start re-ranking'.format(t2 - t1))
        t1 = time.perf_counter()
        psm = fdr(psm)
        t2 = time.perf_counter()
        print('{} secs to finish reranking, ready to write file'.format(t2 - t1))
        '''load protein name'''

        with open("{}.csv".format(spectra_path.split('\\')[-1]), 'w', newline='') as csvfile:
            fieldnames = ['CID_scan', 'mass', 'marker', 're_rank_score', 'alpha', 'pos_a', 'a1_score', 'a2_score',
                          'a_mass', 'a_protein', 'beta', 'pos_b', 'b1_score', 'b2_score', 'b_mass', 'b_protein',
                          'q_value']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for spectrum in psm:
                a_protein = [protein_name[i][0].split(' ')[0] for i in set(spectrum['alpha1'][3])]
                a_protein = ';'.join(a_protein)
                b_protein = [protein_name[i][0].split(' ')[0] for i in set(spectrum['beta1'][3])]
                b_protein = ';'.join(b_protein)
                writer.writerow(
                    {'CID_scan': int(spectrum['CID_scan']), 'mass': spectrum['mass'], 'marker': spectrum['marker'],
                     're_rank_score': spectrum['re_rank'], 'alpha': spectrum['alpha1'][1],
                     'pos_a': spectrum['alpha1'][2],
                     'a1_score': spectrum['alpha1'][0], 'a_mass': spectrum['alpha1'][4], 'a_protein': a_protein,
                     'beta': spectrum['beta1'][1], 'pos_b': spectrum['beta1'][2], 'b1_score': spectrum['beta1'][0],
                     'b_protein': b_protein, 'q_value': spectrum['q_value'], 'b_mass': spectrum['beta1'][4],
                     'a2_score': spectrum['alpha2'], 'b2_score': spectrum['beta2']})
        print('finished file writing {}'.format(spectra_path.split('\\')[-1]))
        LOG_FORMAT = "%(asctime)s=====%(levelname)s++++++%(message)s"
        logging.basicConfig(filename="database.log", level=logging.INFO, format=LOG_FORMAT)
        logging.info("Elapsed {}s.".format(time_diff))
