import os


def decoy_generation(path):
    """
    Generate decoy database (fasta) by revering the target sequence
    :param path:
    """
    isExists = os.path.exists(path.split('.')[0] + '_decoy' + '.fasta')  # check if the decoy already exists
    if isExists:
        os.remove(path.split('.')[0] + '_decoy' + '.fasta')
    _db = read_fasta(path)
    with open(path.split('.')[0] + '_decoy' + '.fasta', 'w') as target_decoy_file:
        for _des, _seq in _db:

            target_decoy_file.write('>' + _des + '\n')
            target_decoy_file.write(''.join([('%s\n' % _seq[i:i + 70])
                                             for i in range(0, len(_seq), 70)]) + '\n')

            target_decoy_file.write('>' + 'DECOY_' + _des + '\n')
            _seq = _seq[::-1]
            target_decoy_file.write(''.join([('%s\n' % _seq[i:i + 70])
                                             for i in range(0, len(_seq), 70)]) + '\n')

    return path.split('.')[0] + '_decoy' + '.fasta'


def read_fasta(path):
    """
    Parse fasta file and yield the tuple of description and sequence
    :param path:
    """
    with open(path, 'r') as fasta_file:
        first_line_flag = True
        line = fasta_file.readline()
        while line:
            if line[0] != '>':
                temp_sequence_list.append(line.strip())
            elif line[0] == '>' and not first_line_flag:
                yield temp_description, ''.join(temp_sequence_list)
                temp_sequence_list = []
                temp_description = line[1:].strip()

            else:
                first_line_flag = False
                temp_sequence_list = []
                temp_description = line[1:].strip()
            line = fasta_file.readline()
        yield temp_description, ''.join(temp_sequence_list)
