from pyteomics import mass, parser
import time


def fragment(peptide_sequence, xl_pos, chain_mass, pre_mass, update_mass_dict, types=('b', 'y')):
    """
    The function generates all possible de charged ions for fragments of types
    `types`, default setting is b/y ions.
    returns a generator
    """
    shift_mass = pre_mass - chain_mass
    gen = [a for a in range(len(peptide_sequence)) if str.isupper(peptide_sequence[a])]
    for idx in range(len(gen)-1):
        if idx < xl_pos:
            for ion_type in types:
                if ion_type in 'abc':
                    yield mass.fast_mass2(
                        peptide_sequence[:gen[idx] + 1], ion_type=ion_type, aa_mass=update_mass_dict)
                else:
                    yield mass.fast_mass2(
                        peptide_sequence[gen[idx] + 1:], ion_type=ion_type, aa_mass=update_mass_dict) + shift_mass
        else:
            for ion_type in types:
                if ion_type in 'abc':
                    yield mass.fast_mass2(
                        peptide_sequence[:gen[idx] + 1], ion_type=ion_type, aa_mass=update_mass_dict) + shift_mass
                else:
                    yield mass.fast_mass2(
                        peptide_sequence[gen[idx] + 1:], ion_type=ion_type, aa_mass=update_mass_dict)


def fast_b_y_fragment(peptide_sequence, xl_pos, chain_mass, pre_mass, update_mass_dict):
    """
        The function generates all possible de charged ions for fragments of b/y ions.
        returns sorted list
        """
    res = []
    shift_mass = pre_mass - chain_mass
    gen = [a for a in range(len(peptide_sequence)) if str.isupper(peptide_sequence[a])]
    first_b = peptide_sequence[0:gen[0] + 1]

    if len(first_b) == 1:
        current_b_ion = update_mass_dict[first_b]
    else:
        current_b_ion = update_mass_dict[first_b[1:-2]] + update_mass_dict[first_b[-1]]

    res.append(current_b_ion) if xl_pos != 0 else res.append(current_b_ion + shift_mass)
    res.append(pre_mass - current_b_ion) if xl_pos != 0 else res.append(chain_mass - current_b_ion)
    for idx in range(1, len(gen)-1):
        parse_b = peptide_sequence[gen[idx-1] + 1: gen[idx] + 1]
        current_b_ion += update_mass_dict[parse_b] if len(parse_b) == 1 else \
            update_mass_dict[parse_b[1:-2]] + update_mass_dict[parse_b[-1]]
        if idx < xl_pos:
            res.append(current_b_ion)
            res.append(pre_mass - current_b_ion)
        else:
            res.append(current_b_ion + shift_mass)
            res.append(chain_mass - current_b_ion)
    res.sort()
    return res


def fast_b_y_c_z_fragment(peptide_sequence, xl_pos, chain_mass, pre_mass, update_mass_dict):
    """
        The function generates all possible de charged ions for fragments of c/z ions.
        returns sorted b/y and c/z list
        """
    by_res = []
    cz_res = []
    c_mass = 17.026549  # mass of NH3
    shift_mass = pre_mass - chain_mass
    gen = [a for a in range(len(peptide_sequence)) if str.isupper(peptide_sequence[a])]
    first_b = peptide_sequence[0:gen[0] + 1]

    if len(first_b) == 1:
        current_b_ion = update_mass_dict[first_b]
        current_c_ion = update_mass_dict[first_b] + c_mass
    else:
        current_b_ion = update_mass_dict[first_b[1:-2]] + update_mass_dict[first_b[-1]]
        current_c_ion = update_mass_dict[first_b[1:-2]] + update_mass_dict[first_b[-1]] + c_mass

    if xl_pos != 0:
        by_res.append(current_b_ion)
        by_res.append(pre_mass - current_b_ion)
        cz_res.append(current_c_ion)
        cz_res.append(pre_mass - current_c_ion)
    else:
        by_res.append(current_b_ion + shift_mass)
        by_res.append(chain_mass - current_b_ion)
        cz_res.append(current_c_ion + shift_mass)
        cz_res.append(chain_mass - current_c_ion)
    for idx in range(1, len(gen)-1):
        parse_b = peptide_sequence[gen[idx-1] + 1: gen[idx] + 1]
        current_b_ion += update_mass_dict[parse_b] if len(parse_b) == 1 else \
            update_mass_dict[parse_b[1:-2]] + update_mass_dict[parse_b[-1]]
        current_c_ion += update_mass_dict[parse_b] if len(parse_b) == 1 else \
            update_mass_dict[parse_b[1:-2]] + update_mass_dict[parse_b[-1]]
        if idx < xl_pos:
            by_res.append(current_b_ion)
            by_res.append(pre_mass - current_b_ion)
            cz_res.append(current_c_ion)
            cz_res.append(pre_mass - current_c_ion)
        else:
            by_res.append(current_b_ion + shift_mass)
            by_res.append(chain_mass - current_b_ion)
            cz_res.append(current_c_ion + shift_mass)
            cz_res.append(chain_mass - current_c_ion)
    by_res.sort()
    cz_res.sort()
    return by_res, cz_res


if __name__ == '__main__':

    print(parser.std_amino_acids)
    print(parser.std_cterm)
    print(parser.std_nterm)
    print(parser.std_labels)
    print(mass.std_aa_mass)

    fix_mod = {'car': [57.021464, ['C']]}
    var_mod = {'34': [34.06, ['ctermK', 'P']]}

    new_mass = {}
    for i, j in fix_mod.items():
        new_mass[i] = j[0]

    for i, j in var_mod.items():
        new_mass[i] = j[0]

    aa_mass = dict(mass.std_aa_mass)
    aa_mass.update(new_mass)
    fix_mod = {key: fix_mod[key][1] for key in fix_mod.keys()}
    var_mod = {key: var_mod[key][1] for key in var_mod.keys()}
    print(fix_mod,var_mod)
    print(aa_mass)
    # print(list(parser.isoforms('KPEPTXISCYDEKEK', variable_mods=var_mod, fixed_mods=fix_mod,max_mods=1)))
    start1 = time.time()
    start2 = time.time()








