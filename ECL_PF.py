import numpy as np
from pyteomics import mass
import concurrent.futures
import os
import json
import pickle
import time
from precursor_discovery import find_precursor_ab_mass, db_to_spectra_extraction
from match import spec_pre_process, cid_p_score, cid_etd_p_score
import csv
from fdr import fdr
from protein_score import protein_score_filtering
import logging


class TopCandidate:
    __slots__ = 'alphaSet', 'betaSet', 'marker'

    def __init__(self, alpha1, beta1, marker):
        self.alphaSet = alpha1
        self.betaSet = beta1
        self.marker = marker


def splitPossibleResult(topAlpha, topBeta, marker, scan, pre_mass):
    """
    Split the top alpha and beta peptides into possible combinations. a2, b2 is the average score here,
    return list of dicts: [{'CID_scan': 123123, 'mass': 799, alpha1 : [-1, 'PEPTIDE', -1, ['0'], 799],
    beta1=[...], alpha2=a2 score, beta2= b2 score},{}]
    :param topAlpha:
    :param topBeta:
    :param marker:
    :param scan:
    :param pre_mass:
    """
    res = []  # create cross-linked peptides list
    flag = []  # avoid redundant assignment
    a2 = []  # used to calculate average scores of top alpha hits
    b2 = []  # used to calculate average scores of top beta hits
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
    a2 = sum(a2) / len(a2)  # average scores of the alpha hits
    b2 = sum(b2) / len(b2)  # average scores of the beta hits
    for _res in res:
        _res['alpha2'] = a2
        _res['beta2'] = b2

    return res


def insertion(candidate, topCandidates, topT):
    """
    Insert the current candidates to the top hits if it's score is high enough.
    Candidate format = (score, sequence, xl position, description, candidate mass),
    topCandidates ordered by the score from smaller to larger, insert candidate to the proper position,
    return the topCandidates with same length if possible
    :param candidate:
    :param topCandidates:
    :param topT:
    """
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
    """
    Searching the spectra in the sub thread
    :param tuple_bag:
    """
    sub_spectra, tol1, tol2, xl, ml, ms, link_site, mass_dict = tuple_bag  # import parameters
    '''find alpha beta chain mass of sub spectra
    returns [{'CID_scan':_, 'mass':_, 'alpha1':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'beta1':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'alpha2':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'beta2':[score, seq, pos, [des1,des2,...], mass, ETD score],
    'marker':_,},{},{}]'''
    t_alpha = 30  # top number of alpha candidates
    t_beta = 30  # top number of beta candidates
    truncation = 0.6  # score requirements, serves like delta score
    res = []  # return candidates
    correspond_matrix = []  # a matrix to find the peptide given their mass in a fast way

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
        process_CID_spec = spec_pre_process(preprocessIn)  # processed CID spectra in SEQUEST way
        process_ETD_spec = []  # to be calculated if exist
        etdConsider = False
        if ele['ETD_peaks']:
            process_ETD_spec = spec_pre_process(ele['ETD_peaks'])
            etdConsider = True

        topAlpha = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_alpha  # initialize a top alpha candidate
        topBeta = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_beta  # initialize a top beta candidate
        marker = 0
        # hit(s) with score,sequence,position,[description],mass

        if etdConsider:  # if CID-ETD concatenated spectra
            for CID_sub_candidate in candidate_list[idx]:  # [beta,alpha,#,int], 3rd layer loop

                if CID_sub_candidate[0] and CID_sub_candidate[1]:
                    specMassRegionTopAlpha = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_alpha
                    specMassRegionTopBeta = [(-1e-3, 'PEPTIDE', -1, '0', 799)] * t_beta

                    for con_alpha in CID_sub_candidate[1]:  # concatenated alpha sequence and description
                        alpha_seq_uncompress = con_alpha[1].split('|')  # alpha peptide sequences
                        alpha_des_uncompress = con_alpha[2].split('|')  # alpha protein names
                        alpha_num = len(alpha_seq_uncompress)  # total possible alpha peptide numbers
                        for idx_con in range(alpha_num):  # 1st layer loop
                            alpha = (con_alpha[0], alpha_seq_uncompress[idx_con], alpha_des_uncompress[idx_con])  # mass, sequence, description

                            '''XlinkX score to match the peptide'''
                            alpha_candidate = cid_etd_p_score(process_CID_spec, process_ETD_spec, alpha[1], ele['mass'],
                                                              alpha[0], tol2, mass_dict, link_site, xl, proN=True)  # return (score, sequence, xl position)
                            if alpha_candidate:
                                currAlpha = (alpha_candidate[0], alpha_candidate[1], alpha_candidate[2], alpha[2], alpha[0])
                                specMassRegionTopAlpha = insertion(currAlpha, specMassRegionTopAlpha, t_alpha)

                    for con_beta in CID_sub_candidate[0]:  # concatenate beta sequence and description, 2nd layer loop
                        beta_seq_uncompress = con_beta[1].split('|')  # beta peptide sequences
                        beta_des_uncompress = con_beta[2].split('|')  # beta protein names
                        beta_num = len(beta_seq_uncompress)  # total possible beta peptide numbers
                        for idx_con in range(beta_num):  # 1st layer loop
                            beta = (con_beta[0], beta_seq_uncompress[idx_con], beta_des_uncompress[idx_con])
                            # mass, sequence, description

                            '''XlinkX score to match the peptide'''
                            beta_candidate = cid_etd_p_score(process_CID_spec, process_ETD_spec, beta[1], ele['mass'],
                                                             beta[0], tol2, mass_dict, link_site, xl, proN=True)  # return (score, sequence, xl position)
                            if beta_candidate:
                                currBeta = (beta_candidate[0], beta_candidate[1], beta_candidate[2], beta[2], beta[0])
                                specMassRegionTopBeta = insertion(currBeta, specMassRegionTopBeta, t_beta)

                    '''compare the score with the highest'''
                    if specMassRegionTopBeta[-1][0] + specMassRegionTopAlpha[-1][0] > topAlpha[-1][0] + topBeta[-1][0]:
                        topAlpha = specMassRegionTopAlpha
                        topBeta = specMassRegionTopBeta
                        marker = CID_sub_candidate[2]

        else:  # only with CID/HCD spectra
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
                            alpha_candidate = cid_p_score(process_CID_spec, alpha[1], ele['mass'], alpha[0], tol2,
                                                          mass_dict, link_site, xl, proN=True)
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
                            beta_candidate = cid_p_score(process_CID_spec, beta[1], ele['mass'], beta[0], tol2,
                                                         mass_dict, link_site, xl, proN=True)
                            # (score, sequence, xl position)
                            if beta_candidate:
                                currBeta = (beta_candidate[0], beta_candidate[1], beta_candidate[2], beta[2], beta[0])
                                specMassRegionTopBeta = insertion(currBeta, specMassRegionTopBeta, t_beta)

                    '''compare the score with the highest'''
                    if specMassRegionTopBeta[-1][0] + specMassRegionTopAlpha[-1][0] > topAlpha[-1][0] + topBeta[-1][0]:
                        topAlpha = specMassRegionTopAlpha
                        topBeta = specMassRegionTopBeta
                        marker = CID_sub_candidate[2]

        '''truncate the top candidates'''
        topAlpha = [_alpha for _alpha in topAlpha if np.isclose(_alpha[0], topAlpha[-1][0], atol=0, rtol=truncation)]
        topBeta = [_beta for _beta in topBeta if np.isclose(_beta[0], topBeta[-1][0], atol=0, rtol=truncation)]
        '''add to the result'''
        res.append((topAlpha, topBeta, marker, ele['CID_scan'], ele['mass']))

    return res


def all_abs_path_of_components(path=r'data\pickled_pair_wise'):
    """
    Retrieve all the spectra file in the file
    :param path:
    """
    res = []
    for maindir, subdir, file_name_list in os.walk(path):
        for filename in file_name_list:
            apath = os.path.join(maindir, filename)
            res.append(apath)

    return res


if __name__ == '__main__':
    """Start the cleavable searching program"""
    with open('ECLPF_conf', 'r') as conf_file:  # load the parameters
        conf = json.load(conf_file)

    Link_site = conf['link_site']  # support multiple sites such as ['K', 'R']
    Fix_mod = conf['fix_mod']  # fixed modification
    Var_mod = conf['var_mod']  # variable modification
    MS1_TOLERANCE = conf['ms1_tol']  # ms1 tolerance
    MS2_TOLERANCE = conf['ms2_tol']  # ms2 tolerance
    xl_mass = conf['xl_mass']  # cross linker residual mass
    m_short = conf['m_short']  # short residual mass
    m_long = conf['m_long']  # longer residual mass

    '''add modification masses into the dictionaries'''
    new_mass = {}  # modifications
    nTerm = {}  # n-term modification dictionary
    cTerm = {}  # c-term modification dictionary

    for i, j in Var_mod.items():  # add variable mods into dictionary
        new_mass[i] = j[0]
        if ('Peptide-nterm' in j[1]) or ('Protein-nterm' in j[1]):
            nTerm[i] = j[0]
        if ('Peptide-cterm' in j[1]) or ('Protein-cterm' in j[1]):
            cTerm[i] = j[0]

    for i, j in Fix_mod.items():  # add fix mods into dictionary
        new_mass[i] = j[0]
        for ele in nTerm:
            new_mass[ele + ')(' + i] = j[0] + nTerm[ele]
        for ele in cTerm:
            new_mass[i + ')(' + ele] = j[0] + nTerm[ele]

    for i, j in Var_mod.items():  # add terminal mods into dictionary
        for ele in nTerm:
            new_mass[ele + ')(' + i] = j[0] + nTerm[ele]
        for ele in cTerm:
            new_mass[i + ')(' + ele] = j[0] + nTerm[ele]

    AA_mass = dict(mass.std_aa_mass)
    AA_mass.update(new_mass)

    '''import spectra file'''
    spectra_list = all_abs_path_of_components()
    print(spectra_list)
    for spectra_path in spectra_list:
        with open(spectra_path, 'rb') as f:
            spectra = pickle.load(f)
        spectra_div = 500  # divided the spectra multiple batches for the multi-processing

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

        with concurrent.futures.ProcessPoolExecutor(max_workers=int(conf['thread'])) as executor:
            result = executor.map(run, args)  # multi-processing

        t2 = time.perf_counter()
        time_diff = round(t2 - t1, 3)
        print(f"It took {time_diff} Secs to execute this method")
        t1 = time.perf_counter()

        '''Start protein feedback processing'''
        with open(r'protein_name', 'rb') as f:
            protein_name = pickle.load(f)
        numberK = {key: protein_name[key][1] for key in protein_name.keys()}
        psm = protein_score_filtering(result, MS1_TOLERANCE, xl_mass, numberK)  # re-ranked CSMs
        t2 = time.perf_counter()
        print('{} seconds to filter the result by protein score, start re-ranking'.format(t2 - t1))
        t1 = time.perf_counter()
        psm = fdr(psm)
        t2 = time.perf_counter()
        print('{} secs to finish reranking, ready to write file'.format(t2 - t1))

        '''Write raw output file without separate FDR control'''
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
        '''Write log file'''
        LOG_FORMAT = "%(asctime)s=====%(levelname)s++++++%(message)s"
        logging.basicConfig(filename="database.log", level=logging.INFO, format=LOG_FORMAT)
        logging.info("Elapsed {}s.".format(time_diff))
