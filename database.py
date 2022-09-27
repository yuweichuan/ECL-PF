import datetime
import time
import pickle
import concurrent.futures
import numpy as np
import os
from decoy_generation import decoy_generation, read_fasta
import shutil
import json
import logging
from collections import deque
import itertools as it
import re

std_aa_mass = {
    'G': 57.02146,
    'A': 71.03711,
    'S': 87.03203,
    'P': 97.05276,
    'V': 99.06841,
    'T': 101.04768,
    'C': 103.00919,
    'L': 113.08406,
    'I': 113.08406,
    'N': 114.04293,
    'D': 115.02694,
    'Q': 128.05858,
    'K': 128.09496,
    'E': 129.04259,
    'M': 131.04049,
    'H': 137.05891,
    'F': 147.06841,
    'U': 150.95364,
    'R': 156.10111,
    'Y': 163.06333,
    'W': 186.07931,
    'O': 237.14773,
}
expasy_rules = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3': r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                    r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                    r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV])',
    'thrombin': r'((?<=G)R(?=G))|'
                r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
}
_nterm_mod = r'[^-]+-$'
_cterm_mod = r'-[^-]+$'
std_nterm = 'H-'
std_cterm = '-OH'


def mass_calculation(sequence, aa_mass):
    """calculate the mono-isotopic mass given the specific peptide"""
    gen = [a for a in range(len(sequence)) if str.isupper(sequence[a])]
    if sequence[0:gen[0] + 1] not in aa_mass:  # mass of water and the first amino acid
        _mass = 18.01055 + aa_mass[sequence[0:gen[0]]] + aa_mass[sequence[gen[0]]]
    else:
        _mass = 18.01055 + aa_mass[sequence[gen[0]]]
    for _index in range(1, len(gen)):
        if sequence[gen[_index - 1] + 1:gen[_index] + 1] in aa_mass:
            _mass += aa_mass[sequence[gen[_index]]]
        else:
            _mass += aa_mass[sequence[gen[_index - 1] + 1:gen[_index]]] + aa_mass[sequence[gen[_index]]]
    return _mass


def cleave(sequence, rule, missed_cleavages=0, min_length=1, Methionine_drop=True):
    """cleave the protein into peptides and specify their relative positions"""
    rule = expasy_rules.get(rule, rule)
    peptides = []
    ml = missed_cleavages + 2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1

    for ii in it.chain([x.end() for x in re.finditer(rule, sequence)], [None]):
        cleavage_sites.append(ii)
        if cl < ml:
            cl += 1
        for jj in trange[:cl - 1]:
            seq = sequence[cleavage_sites[jj]:cleavage_sites[-1]]
            if seq and len(seq) >= min_length:
                peptides.append((seq, cleavage_sites[jj]))  # peptide sequence and position (index from 0)
    if Methionine_drop:
        '''add methionine dropped peptides'''
        meth_group = [ele[1:] for ele, index in peptides if (ele[0] == 'M' and index == 0 and len(ele) > min_length)]
        for _ele in meth_group:
            peptides.append((_ele, 1))
    return peptides


def is_term_mod(label):
    return (re.match(_nterm_mod, label) or re.match(_cterm_mod, label)) is not None


def pep_parse(sequence):
    sequence = str(sequence)
    if len(sequence) >= 2:
        return [(std_nterm, sequence[0])] + [tuple(ele) for ele in sequence[1:-1]] + [(sequence[-1], std_cterm)]
    elif len(sequence) == 1:
        return [(std_nterm, sequence, std_cterm)]
    else:
        raise ValueError("not a valid peptide sequence")


def tostring(parsed_sequence):
    """return sequence and dict of additional num of PTMs"""
    parsed_sequence = list(parsed_sequence)
    labels = []
    nterm = parsed_sequence[0]
    cterm = parsed_sequence[-1]
    ptmNum = dict()

    for ele in (m for m in parsed_sequence[1:-1] if len(m) > 1):
        ptmNum[ele[0][1: -1]] = ptmNum.setdefault(ele[0][1: -1], 0) + 1

    labels.extend(''.join(g) for g in parsed_sequence[1: -1])
    return ''.join(labels), ptmNum


def isoforms(sequence, **kwargs):
    """
    Apply variable and fixed modifications to the polypeptide and yield
    the unique modified sequences.
    out : iterator over strings or lists
        All possible unique polypeptide sequences resulting from
        the specified modifications are yielded one by one.
    """

    def main(_group):  # index of the residue (capital letter) in `group`
        if _group[-1][0] == '-':
            ind = -2
        else:
            ind = -1
        return len(_group) + ind, _group[ind]

    def apply_mod(label, mod):
        _group = list(label)
        m = main(_group)[0]
        if m == 0 and not is_term_mod(mod):
            _group.insert(0, '(' + mod + ')')
        elif mod[0] == '-' and (_group[-1] == std_cterm):
            _group[-1] = mod
        elif mod[-1] == '-' and (_group[0] == std_nterm):
            _group[0] = mod
        elif not is_term_mod(mod):
            if m and not _group[m - 1][-1] == '-':
                pass
            else:
                _group.insert(m, '(' + mod + ')')
        return tuple(_group)

    variable_mods = kwargs.get('variable_mods', {})
    nTermMark = kwargs.get('ProteinN', False)
    cTermMark = kwargs.get('ProteinC', False)
    peptideN = []
    peptideC = []
    proteinN = []
    proteinC = []
    for k, v in variable_mods.items():
        if 'Protein-nterm' in v:
            proteinN.append(k)

        if 'Protein-cterm' in v:
            proteinC.append(k)

        if 'Peptide-nterm' in v:
            peptideN.append(k)

        if 'Peptide-cterm' in v:
            peptideC.append(k)

    fixed_mods = kwargs.get('fixed_mods', {})
    parsed = pep_parse('#' + sequence + '#')  # append the sequence
    max_mods = kwargs.get('max_mods')

    # Apply fixed modifications
    for cmod in fixed_mods:
        for ii, group in enumerate(parsed):
            if main(group)[1] in fixed_mods[cmod]:
                parsed[ii] = apply_mod(group, cmod)

    # Create a list of possible states for each group
    # Start with N-terminal mods and regular mods on the N-terminal residue
    second = set(apply_mod(parsed[0], m) for m, r in variable_mods.items()
                 if (r is True or
                     main(parsed[0])[1] in r or
                     'nterm' + main(parsed[0])[1] in r or
                     (len(parsed) == 1 and 'cterm' + main(parsed[0])[1] in r))
                 and not is_term_mod(m)
                 ).union([parsed[0]])
    first = it.chain((apply_mod(group, mod) for group in second
                      for mod, res in variable_mods.items()
                      if (mod.endswith('-') or (mod.startswith('-') and len(parsed) == 1))
                      and (res is True or main(group)[1] in res)), second)
    states = [[parsed[0]] + list(set(first).difference({parsed[0]}))]
    # Continue with regular mods
    states.extend([group] + list(set(apply_mod(group, mod)
                                     for mod in variable_mods if (
                                             variable_mods[mod] is True or
                                             group[-1] in variable_mods[mod]) and not is_term_mod(mod)
                                     ).difference({group}))
                  for group in parsed[1:-1])
    # Finally add C-terminal mods and regular mods on the C-terminal residue
    if len(parsed) > 1:
        second = set(apply_mod(parsed[-1], m) for m, r in variable_mods.items()
                     if (r is True or
                         main(parsed[-1])[1] in r or
                         'cterm' + main(parsed[-1])[1] in r)
                     and not is_term_mod(m)
                     ).union((parsed[-1],))
        first = it.chain((apply_mod(group, mod) for group in second
                          for mod, res in variable_mods.items()
                          if mod.startswith('-') and (
                                  res is True or main(group)[1] in res)), second)
        states.append([parsed[-1]] + list(set(first).difference({parsed[-1]})))
    sites = [s for s in enumerate(states) if len(s[1]) > 1]
    if max_mods is None or max_mods > len(sites):
        possible_states = it.product(*states)
    else:
        def state_lists():
            for m in range(max_mods + 1):
                for comb in it.combinations(sites, m):
                    skel = [[s[0]] for s in states]
                    for ind, e in comb:
                        skel[ind] = e[1:]
                    yield skel
        possible_states = it.chain.from_iterable(
            it.product(*skel) for skel in state_lists())
    possible_states = [i for i in possible_states]

    res = [tostring(form) for form in possible_states]
    updateProN = []
    updateProC = []
    updatePepN = []
    updatePepC = []
    if nTermMark and proteinN:
        for ele in proteinN:
            for seqRes, dictRes in res:
                newSeq = '(' + ele + ')' + str(seqRes)
                newDict = dict(dictRes)
                newDict[ele] = newDict.setdefault(ele, 0) + 1
                newEle = tuple((newSeq, newDict))
                updateProN.append(newEle)
    if cTermMark and proteinC:
        for ele in proteinC:
            for seqRes, dictRes in res:
                newSeq = str(seqRes)[:-1] + '(' + ele + ')' + str(seqRes)[-1]
                newDict = dict(dictRes)
                newDict[ele] = newDict.setdefault(ele, 0) + 1
                newEle = tuple((newSeq, newDict))
                updateProC.append(newEle)

    for ele in peptideN:
        for seqRes, dictRes in res:
            newSeq = '(' + ele + ')' + str(seqRes)
            newDict = dict(dictRes)
            newDict[ele] = newDict.setdefault(ele, 0) + 1
            newEle = tuple((newSeq, newDict))
            updatePepN.append(newEle)

    for ele in peptideC:
        for seqRes, dictRes in res:
            newSeq = str(seqRes)[:-1] + '(' + ele + ')' + str(seqRes)[-1]
            newDict = dict(dictRes)
            newDict[ele] = newDict.setdefault(ele, 0) + 1
            newEle = tuple((newSeq, newDict))
            updatePepC.append(newEle)
    return res + updateProN + updateProC + updatePepN + updatePepC


def database_generate(tuple_bag, interval_dalton=100, maxda=6000):
    """primary amine, peptide n/c and protein n/c needed"""
    fas_list, parse_rule, num_max_mod, max_length, min_length, miss_cleavage,\
        link_site, var_mod, fix_mod, aa_mass, num, div = tuple_bag
    db_peptides = {}
    number = int(num) * div
    for description, sequence in fas_list:
        number += 1
        proteinLength = len(sequence)

        if all(word not in ['B', 'J', 'X', 'Z', 'O', 'U'] for word in sequence):
            # form_description = description.split(' ')[0]
            form_description = number
            new_peptides = cleave(sequence, parse_rule, missed_cleavages=miss_cleavage, min_length=min_length)

            # new_peptides = set([ele[0] for ele in new_peptides])  # temporarily
            new_peptides = (new_peptide for new_peptide in new_peptides
                            if len(new_peptide[0]) <= max_length and
                            (not set(new_peptide[0][0:-1]).isdisjoint(link_site) or new_peptide[1] <= 1))
            # include the n-term site and primary amine by default
            for element, index in new_peptides:
                if index <= 1:
                    proteinN = True
                else:
                    proteinN = False
                if len(element) + index == proteinLength:
                    proteinC = True
                else:
                    proteinC = False
                forms = isoforms(element, variable_mods=var_mod, fixed_mods=fix_mod, max_mods=num_max_mod,
                                 ProteinN=proteinN, ProteinC=proteinC)
                # forms is a generator
                baseMass = mass_calculation(element, aa_mass=aa_mass)

                for form, ptmDict in forms:
                    form_mass = round(baseMass + sum([aa_mass[ii] * jj for ii, jj in ptmDict.items()]), 3)
                    # precision 3 decimal is enough!!!

                    # if link_site in (form)[0:-1] and len(parser.parse(form)[-1]) == 1 and form_mass < maxda:
                    # include the n-term site, modified site cannot link and last site cannot mod

                    # if link_site in parser.parse(form)[0:-1] and form_mass < maxda:
                    # include the n-term site, modified site cannot link

                    if form_mass < maxda:
                        # include the n-term site and c-term site and primary amine by default
                        form = '[' + form if proteinN else form
                        form = form + ']' if proteinC else form
                        if form_mass not in db_peptides.keys():
                            db_peptides[form_mass] = [(form, form_description)]
                        else:
                            db_peptides[form_mass].append((form, form_description))

    # with open('database_file/database_{}'.format(num), 'wb') as pickout:
    #     pickle.dump(db_peptides, pickout)
    for Da in range(0, maxda, interval_dalton):
        sub_db_keys = np.around(np.arange(Da - 0.2, Da + interval_dalton + 0.2, 0.001), 3)
        # 0.2 Da leeway in case of extreme mass
        sub_db = {key: db_peptides[key] for key in db_peptides.keys() & sub_db_keys}

        with open(r'database_file/database_{}_{}_Da'.format(num, Da), 'wb') as pickout:
            pickle.dump(sub_db, pickout)
    '''write log'''
    LOG_FORMAT = "%(asctime)s=====%(levelname)s++++++%(message)s"
    logging.basicConfig(filename="database.log", level=logging.INFO, format=LOG_FORMAT)
    logging.info("Finished the {} th batch data parsing.".format(num))
    del db_peptides


def merge_fasta():
    """merge the same mass interval database together and provide a mass index set
    returns: [set(mass1,mass2,...), dict(mass1:[(form,des),(form,des),...], mass2:[(form,des),(form,des),..])]"""
    os.chdir('database_file')
    files = os.listdir()
    file_key = {name.split('_')[-2] for name in files}
    for k in file_key:
        file = [name for name in files if name.split('_')[-2] == k]
        merge_list = [set(), dict()]
        for fi in file:
            with open(fi, 'rb') as f:
                sub_merge = pickle.load(f)
            for key in sub_merge.keys():
                if key in merge_list[1].keys():
                    merge_list[1][key] = merge_list[1][key] + sub_merge[key]
                else:
                    merge_list[1][key] = sub_merge[key]
                    merge_list[0].add(key)
        tempDict = dict()
        for _key, _value in merge_list[1].items():
            tempValue = dict()
            tempMarker = set()
            for ele_value in _value:
                if ele_value[0] not in tempMarker:
                    tempMarker.add(ele_value[0])
                    tempValue[ele_value[0]] = str(ele_value[1])
                else:
                    tempValue[ele_value[0]] = tempValue[ele_value[0]] + '$' + str(ele_value[1])
            concatenateSequence = '|'.join(tempValue.keys())
            concatenateDescription = '|'.join(tempValue.values())
            tempDict[_key] = (concatenateSequence, concatenateDescription)
        merge_list[1] = tempDict

        with open('{}'.format(k), 'wb') as f:
            pickle.dump(merge_list, f)
    for path in files:
        os.remove(path)


if __name__ == '__main__':
    '''import configuration'''
    with open('ECLPF_conf', 'r') as conf_file:
        conf = json.load(conf_file)
    Parse_rule = conf['parse_rule']
    Fasta_path = conf['fasta_path']
    Max_length = conf['max_length']
    Min_length = conf['min_length']
    Miss_cleavage = conf['miss_cleavage']
    Num_max_mod = conf['num_max_mod']  # exclude fixed modifications
    Link_site = conf['link_site']  # support multiple sites such as ['K', 'R']
    Fix_mod = conf['fix_mod']
    Var_mod = conf['var_mod']

    database_file_isExists = os.path.exists('database_file')
    if database_file_isExists:
        shutil.rmtree('database_file')
        time.sleep(2)
        os.mkdir('database_file')
    else:
        os.mkdir('database_file')

    starttime = datetime.datetime.now()
    '''add modification masses into the dict'''
    new_mass = {}
    for i, j in Fix_mod.items():
        new_mass[i] = j[0]

    for i, j in Var_mod.items():
        new_mass[i] = j[0]

    AA_mass = std_aa_mass
    AA_mass.update(new_mass)
    Fasta_list = list(read_fasta(decoy_generation(Fasta_path)))  # adding decoy database
    # Fasta_list = list(fasta.read(Fasta_path))  # without decoy database
    '''construct the mapping of protein name and index'''
    protein_name = dict()
    upperScore = 100
    protein_name[0] = ('PEPTIDE', upperScore - 1)
    name_index = 0
    for des, seq in Fasta_list:
        name_index += 1
        kk = min(upperScore - 1, seq.count('K') + 1)
        protein_name[name_index] = (des, kk)
    with open(r'protein_name', 'wb') as f:
        pickle.dump(protein_name, f)

    '''create the new fix_mod and var_mod mapping'''
    Fix_mod = {key: Fix_mod[key][1] for key in Fix_mod.keys()}
    Var_mod = {key: Var_mod[key][1] for key in Var_mod.keys()}

    div_data = 500  # divide fasta by every # number
    section = [i for i in range(len(Fasta_list)) if i % div_data == 0]
    section_length = len(section)
    num_magnitude = len(str(section_length))
    args = []
    for i in range(section_length):
        args.append((Fasta_list[section[i]:section[i] + div_data], Parse_rule, Num_max_mod, Max_length, Min_length,
                     Miss_cleavage, Link_site, Var_mod, Fix_mod, AA_mass, str(i).zfill(num_magnitude), div_data))
    print('start multiprocessing')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(database_generate, args)

    endtime = datetime.datetime.now()
    print('{} seconds!'.format((endtime - starttime).seconds))
    print('merging data...')
    a1 = time.perf_counter()
    merge_fasta()
    a2 = time.perf_counter()
    print('{} seconds! finished!'.format(a2-a1))
