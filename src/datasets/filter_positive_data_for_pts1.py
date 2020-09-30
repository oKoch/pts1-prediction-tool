from Bio import SeqIO


def read_in_aminoacids(rel_path, aacid_length, is_targeted):
    raw = list(SeqIO.parse(rel_path, "fasta"))
    prepared_seq = []
    targeted = []
    for entry in raw:
        aminoacid_seq = list(str(entry.seq[aacid_length:]).upper())
        targeted.append(is_targeted)
        prepared_seq.append(aminoacid_seq)

    return prepared_seq, targeted


minus_3 = ['S', 'A', 'C', 'P', 'H', 'T', 'N', 'Q', 'E', 'G', 'V']
minus_2 = ['K', 'R', 'H', 'Q', 'D', 'N', 'S', 'M']
minus_1 = ['L', 'F', 'I', 'M', 'Y']


def read_in_aminoacids_with_char_to_remove(rel_path):
    raw = list(SeqIO.parse(rel_path, "fasta"))
    aacid_length = -3
    prepared_seq = []
    entries_with_potential_pts1 = []
    for entry in raw:
        is_ok = False
        seq = str(entry.seq).upper()
        aminoacid_seq = list(seq.replace('*', '')[aacid_length:])
        if str(aminoacid_seq[-3]) in minus_3:
            if str(aminoacid_seq[-2]) in minus_2:
                if str(aminoacid_seq[-1]) in minus_1:
                    is_ok = True
        if is_ok:
            entries_with_potential_pts1.append(entry)

    return entries_with_potential_pts1


min_length_aminoacids = -15
dicty_positiv_seq = read_in_aminoacids_with_char_to_remove("./raw_peroxisomal_sequences.fasta")


with open("pts1_filtered_peroxisomal_sequences.fasta", "w") as output_handle:
    SeqIO.write(dicty_positiv_seq, output_handle, "fasta")

