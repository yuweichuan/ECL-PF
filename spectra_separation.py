from pyteomics import mzxml
import os
import shutil
import pickle
import concurrent.futures
import time
import numpy as np
import json


def separate(tuple_bag):
    """
    Separate the spectra by CID and ETD type. And merge the peaks within the MS tolerance
    :param tuple_bag:
    """
    file_name, output, activation_type, tol1, tol2, ratio = tuple_bag  # tol unit is ppm
    with mzxml.read(file_name) as spectra:
        CID = set()
        ETD = set()
        for spectrum in spectra:
            if spectrum['msLevel'] == 2:
                if 'activationMethod' not in spectrum['precursorMz'][0].keys():
                    activation = activation_type[0]
                else:
                    activation = spectrum['precursorMz'][0]['activationMethod']

                if activation == activation_type[0]:
                    CID.add(spectrum['num'])
                elif activation == activation_type[1]:
                    ETD.add(spectrum['num'])

    with open(output, 'r') as lines:  # output is the hardklor deconvolution (hk) file path
        CID_spectra = dict()
        ETD_spectra = dict()  # scan:{'mass': float, 'peaks':[(mass1, intensity1),(mass2, intensity2),...]}
        CIDETD_spectrum = []
        for line in lines:
            new_line = line.split()
            if new_line[0] == 'S':
                if CIDETD_spectrum:
                    if str(CIDETD_spectrum[1]) in CID:
                        CID_spectra[CIDETD_spectrum[1]] = {'mz': CIDETD_spectrum[3], 'mass': CIDETD_spectrum[0],
                                                           'charge': CIDETD_spectrum[4], 'peaks': CIDETD_spectrum[2]}
                        # {scan:{'mz':m/z, 'mass':mass, 'peaks':[(peak1, intensity1),(peak2, intensity2)]}}
                    elif str(CIDETD_spectrum[1]) in ETD:
                        ETD_spectra[CIDETD_spectrum[1]] = {'mz': CIDETD_spectrum[3], 'mass': CIDETD_spectrum[0],
                                                           'charge': CIDETD_spectrum[4], 'peaks': CIDETD_spectrum[2]}
                CIDETD_spectrum = [round(float(new_line[4]), 3), int(new_line[1]), [], round(float(new_line[6]), 3),
                                   int(new_line[5])]
                # make sure mass decimal is 4
                # [precursor mass, scan, empty peak list (mass, intensity), mz]
            else:
                CIDETD_spectrum[2].append((round(float(new_line[1]), 3), int(new_line[3])))
                # make sure mass decimal is 4
                # expand peak list with (mass, intensity) tuple
        if str(CIDETD_spectrum[1]) in CID:
            CID_spectra[CIDETD_spectrum[1]] = {'mz': CIDETD_spectrum[3], 'charge': CIDETD_spectrum[4],
                                               'mass': CIDETD_spectrum[0], 'peaks': CIDETD_spectrum[2]}
            # {scan:{'mz':m/z, 'mass':mass, 'peaks':[(peak1, intensity1),(peak2, intensity2)]}}
        elif str(CIDETD_spectrum[1]) in ETD:
            ETD_spectra[CIDETD_spectrum[1]] = {'mz': CIDETD_spectrum[3], 'charge': CIDETD_spectrum[4],
                                               'mass': CIDETD_spectrum[0], 'peaks': CIDETD_spectrum[2]}

    pair_wise = []  # [{'mass':xxx, 'charge': xxx, 'CID_scan':xxx, 'CID_peaks':[(mass1,intensity1),(mass2,intensity2)],
    # 'ETD_peaks:[(mass1,intensity1),(mass2,intensity2)]}]
    for cid in CID_spectra.keys():
        etd_candidates_scan = [etd_scan for etd_scan in range(cid - 10, cid + 10) if str(etd_scan) in ETD]
        etd_candidates = [ETD_spectra[scan] for scan in etd_candidates_scan if
                          np.isclose(ETD_spectra[scan]['mass'], CID_spectra[cid]['mass'], rtol=tol1) and
                          np.isclose(ETD_spectra[scan]['mz'], CID_spectra[cid]['mz'], rtol=tol1)]

        if len(etd_candidates) == 1 and \
                len(CID_spectra[cid]['peaks']) >= max(4, int(CID_spectra[cid]['mass'] / 111) * 2 * ratio) and \
                len(etd_candidates[0]['peaks']) > max(1, int(CID_spectra[cid]['mass'] / 111) * 2 * ratio):
            # add paired CID ETD spectrum
            # the average amino acid mass is 111 Da, require a minimum 4 peaks in each CID spectrum because of 4 markers
            # but not require ETD spectrum with at least 4 peaks

            '''merge the mass with different charges first'''
            CID_peaks = sorted(CID_spectra[cid]['peaks'])
            merged_CID = [CID_peaks[0]]
            for cid_idx in range(1, len(CID_peaks)):
                if np.isclose(merged_CID[-1][0], CID_peaks[cid_idx][0], rtol=tol2):
                    weighted_mass = (merged_CID[-1][0] * merged_CID[-1][1] + CID_peaks[cid_idx][0] *
                                     CID_peaks[cid_idx][1]) / (merged_CID[-1][1] + CID_peaks[cid_idx][1])
                    sum_intensity = merged_CID[-1][1] + CID_peaks[cid_idx][1]
                    merged_CID[-1] = (round(weighted_mass, 3), sum_intensity)
                else:
                    merged_CID.append((CID_peaks[cid_idx][0], CID_peaks[cid_idx][1]))
            ETD_peaks = sorted(etd_candidates[0]['peaks'])
            merged_ETD = [ETD_peaks[0]]
            for etd_idx in range(1, len(ETD_peaks)):
                if np.isclose(merged_ETD[-1][0], ETD_peaks[etd_idx][0], rtol=tol2):
                    weighted_mass = (merged_ETD[-1][0] * merged_ETD[-1][1] + ETD_peaks[etd_idx][0] *
                                     ETD_peaks[etd_idx][1]) / (merged_ETD[-1][1] + ETD_peaks[etd_idx][1])
                    sum_intensity = merged_ETD[-1][1] + ETD_peaks[etd_idx][1]
                    merged_ETD[-1] = (round(weighted_mass, 3), sum_intensity)
                else:
                    merged_ETD.append((ETD_peaks[etd_idx][0], ETD_peaks[etd_idx][1]))
            pair_wise.append({'mass': CID_spectra[cid]['mass'], 'CID_scan': cid, 'charge': CID_spectra[cid]['charge'],
                              'CID_peaks': merged_CID, 'ETD_peaks': merged_ETD})  # sort peaks by their mass

        elif len(etd_candidates) == 0 and \
                len(CID_spectra[cid]['peaks']) >= max(4, int(CID_spectra[cid]['mass'] / 111) * 2 * ratio):
            # add single CID spectrum
            # the average amino acid mass is 111 Da, require a minimum 4 peaks in each CID spectrum because of 4 markers
            CID_peaks = sorted(CID_spectra[cid]['peaks'])
            merged_CID = [CID_peaks[0]]
            for cid_idx in range(1, len(CID_peaks)):
                if np.isclose(merged_CID[-1][0], CID_peaks[cid_idx][0], rtol=tol2):
                    weighted_mass = (merged_CID[-1][0] * merged_CID[-1][1] + CID_peaks[cid_idx][0] *
                                     CID_peaks[cid_idx][1]) / (merged_CID[-1][1] + CID_peaks[cid_idx][1])
                    sum_intensity = merged_CID[-1][1] + CID_peaks[cid_idx][1]
                    merged_CID[-1] = (round(weighted_mass, 3), sum_intensity)
                else:
                    merged_CID.append((CID_peaks[cid_idx][0], CID_peaks[cid_idx][1]))
            pair_wise.append({'mass': CID_spectra[cid]['mass'], 'CID_scan': cid, 'charge': CID_spectra[cid]['charge'],
                              'CID_peaks': merged_CID, 'ETD_peaks': None})  # sort peaks by their mass
    with open('pickled_pair_wise/{}'.format(file_name.split('/')[-1].split('.')[0]), 'wb') as pickout:
        # truncate string
        pickle.dump(pair_wise, pickout)
    del pair_wise


if __name__ == '__main__':
    with open('ECLPF_conf', 'r') as conf_file:
        conf = json.load(conf_file)
    TOLERANCE1 = conf['ms1_tol']  # this is for the use of merge ETD, CID spectra
    TOLERANCE2 = conf['ms2_tol']  # this is for the use of same fragment with different charges, we merge them together
    PEAK_RATIO = 0  # require peak ratio with at least proportion for CID and ETD spectra (domain [0,1])
    os.chdir('data/')
    Activation_Type = conf['activation_type']
    for ele in conf['data_path']:
        if ' ' in ele:
            print('space contains')
            quit()
    File_name = [name.split('/')[-1] for name in conf['data_path']]
    print(File_name)

    Output = [r'hardklor_result/{}'.format(name.replace('mzXML', 'txt')) for name in File_name]
    print(Output)
    hardklor_isExists = os.path.exists('hardklor_result')
    pairwise_isExists = os.path.exists('pickled_pair_wise')
    if pairwise_isExists:
        shutil.rmtree('pickled_pair_wise')
        time.sleep(2)
        os.mkdir('pickled_pair_wise')
    else:
        os.mkdir('pickled_pair_wise')

    if hardklor_isExists:
        shutil.rmtree('hardklor_result')
        time.sleep(2)
        os.mkdir('hardklor_result')
    else:
        os.mkdir('hardklor_result')

    '''Use Hardklor to de-charge and generate txt file'''
    with open('Hardklor.conf') as f:
        param = f.read()

    with open('parameter.conf', 'w') as f:
        f.write(param)
        for i in range(len(conf['data_path'])):
            f.write("\n{} {}".format(conf['data_path'][i].replace('\\', '/'), Output[i]))

    os.system("Hardklor.exe %s" % 'parameter.conf')

    os.remove('parameter.conf')
    args = []

    print('Converting results...')
    for i in range(len(conf['data_path'])):
        args.append((conf['data_path'][i], Output[i], Activation_Type, TOLERANCE1, TOLERANCE2, PEAK_RATIO))

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        a = executor.map(separate, args)

    print('Done!')
