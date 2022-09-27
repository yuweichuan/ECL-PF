from pyteomics import fasta
import os


def decoy_generation(path):  # by reversing the fasta
    isExists = os.path.exists(path.split('.')[0] + '_decoy' + '.fasta')
    if isExists:
        os.remove(path.split('.')[0] + '_decoy' + '.fasta')
    # fasta.write_decoy_db(path, mode='reverse', output=path.split('.')[0] + '_decoy' + '.fasta')
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
    """parse fasta file and yield the tuple of description and sequence"""
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


if __name__ == '__main__':
    a = decoy_generation('BSA_test.fasta')

    sp = list(read_fasta(a))
    for i, j in sp:
        print(i)
        print(j)
