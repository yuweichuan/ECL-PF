def fast_b_y_fragment(peptide_sequence, xl_pos, chain_mass, pre_mass, update_mass_dict):
    """
    The function generates all possible de charged ions
    for fragments of b/y ions. Returns sorted list
    :param peptide_sequence:
    :param xl_pos:
    :param chain_mass:
    :param pre_mass:
    :param update_mass_dict:
    """
    res = []  # Returned list
    shift_mass = pre_mass - chain_mass  # peak shift mass due to the cross-linked counter part.
    gen = [a for a in range(len(peptide_sequence)) if str.isupper(peptide_sequence[a])]  # amino acid list
    first_b = peptide_sequence[0:gen[0] + 1]  # start (b) ions

    if len(first_b) == 1:
        current_b_ion = update_mass_dict[first_b]  # current b series ions
    else:
        current_b_ion = update_mass_dict[first_b[1:-2]] + update_mass_dict[first_b[-1]]

    res.append(current_b_ion) if xl_pos != 0 else res.append(current_b_ion + shift_mass)
    res.append(pre_mass - current_b_ion) if xl_pos != 0 else res.append(chain_mass - current_b_ion)
    for idx in range(1, len(gen)-1):
        parse_b = peptide_sequence[gen[idx-1] + 1: gen[idx] + 1]
        current_b_ion += update_mass_dict[parse_b] if len(parse_b) == 1 else \
            update_mass_dict[parse_b[1:-2]] + update_mass_dict[parse_b[-1]]
        if idx < xl_pos:
            res.append(current_b_ion)  # add b ions
            res.append(pre_mass - current_b_ion)  # add mass shifted y ions
        else:
            res.append(current_b_ion + shift_mass)  # add mass shifted b ions
            res.append(chain_mass - current_b_ion)  # add y ions
    res.sort()
    return res


def fast_b_y_c_z_fragment(peptide_sequence, xl_pos, chain_mass, pre_mass, update_mass_dict):
    """
    The function generates all possible de charged ions for fragments of b/y c/z ions.
    Returns sorted b/y and c/z list.
    :param peptide_sequence:
    :param xl_pos:
    :param chain_mass:
    :param pre_mass:
    :param update_mass_dict:
    """
    by_res = []  # b/y ion list
    cz_res = []  # c/z ion list
    c_mass = 17.026549  # mass of NH3
    shift_mass = pre_mass - chain_mass
    gen = [a for a in range(len(peptide_sequence)) if str.isupper(peptide_sequence[a])]  # amino acid list
    first_b = peptide_sequence[0:gen[0] + 1]  # start (b) ions

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
            by_res.append(current_b_ion)  # add b ions
            by_res.append(pre_mass - current_b_ion)  # add mass shifted y ions
            cz_res.append(current_c_ion)  # add c ions
            cz_res.append(pre_mass - current_c_ion)  # add mass shifted z ions
        else:
            by_res.append(current_b_ion + shift_mass)  # add mass shifted b ions
            by_res.append(chain_mass - current_b_ion)  # add y ions
            cz_res.append(current_c_ion + shift_mass)  # add mass shifted c ions
            cz_res.append(chain_mass - current_c_ion)  # add z ions
    by_res.sort()
    cz_res.sort()
    return by_res, cz_res
