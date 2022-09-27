import json

'''this is the input configuration file'''

conf = {'parse_rule': 'trypsin',  # 'arg-c', 'asp-n', 'pepsin ph1.3', 'pepsin ph2.0' etc.
        'fasta_path': "BSA.fasta",  # path to your fasta file
        'activation_type': ['CID', 'ETD'],  # ['HCD', 'ETD'] or ['CID', 'ETD']
        'max_length': 50,  # maximum peptide length
        'min_length': 5,  # minimum peptide length
        'ms1_tol': 5e-6,  # MS1 mass tolerance in ppm
        'ms2_tol': 2e-5,  # MS2 mass tolerance in ppm
        'miss_cleavage': 2,  # allowed missed cleavage in peptides
        'num_max_mod': 3,  # maximun number of modification allowed, exclude fixed modifications
        'link_site': ['K', '['],  # support multiple sites such as ['K', 'R','[']. '[' means protein N-term
        'fix_mod': {'car': [57.021464, ['C']]},  # fixed modification. 'car' is the modification name, can only use lowercase letter and num

        'var_mod': {'28': [28.0313, ['K', 'Peptide-nterm']],
                    '34': [34.0631, ['K', 'Peptide-nterm']],
                    'oxi': [15.9949, ['M']]},  # variable modification. 'oxi' is the modification name, can only use lowercase letter and num

        'xl_mass': 509.097,  # dsbu 196.0848, dsso 158.0038, dsbso 308.0388, cbdps 509.097
        'm_short': 54.011,  # dsbu 85.0528, dsso 54.0106, dsbso 54.0106, cbdps 54.011
        'm_long': 455.086}  # dsbu 111.032, dsso 85.9826/103.9932, dsbso 236.0177/254.0283, cbdps 455.086 

with open('ECLPF_conf', 'w') as f:
    json.dump(conf, f, indent=4, separators=(',', ':'))
