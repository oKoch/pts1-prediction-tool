import numpy as np
import pandas as pd
import datetime

from pandas.api.types import CategoricalDtype
from Bio import SeqIO, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.svm import SVC

print("Start pts1_prediction_tool")

aminoacid_onelettercode = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                           'W', 'Y', 'X', 'U']
peroxisomal = 1
not_peroxisomal = 0

def prepare_sequence_for_machine_learning(raw, aacid_length, is_targeted):
    aminoacid_sequences = []
    prepared_seq = []
    targeted = []

    for entry in raw:
        # Vereinheitlichung der Aminosäuresequenz zu Großbuchstaben
        seq = str(entry.seq).upper()
        # Aussortieren zu kurzer Aminosäuresequenzen
        if (len(seq) < (-1 * aacid_length)):
            continue

        # Aussortieren der doppelten Aminosäuresequenzen
        if seq not in aminoacid_sequences:
            # Vorbereitung Aussortierung doppelter Sequenzen
            aminoacid_sequences.append(seq)
            # Verkürzen der Sequenz auf letzten stellen vom c-terminus, und umwandlung zu liste, Entfernen des c-terminus sternchen im string
            aminoacid_seq = list(seq.replace('*', '')[aacid_length:])
            targeted.append(is_targeted)
            prepared_seq.append(aminoacid_seq)

    return prepared_seq, targeted


def prepare_dicty_sequences(raw, aacid_length):
    prepared_seq = []
    for entry in raw:
        seq = str(entry.seq).upper()
        if (len(seq) < (-1 * aacid_length)):
            continue

        aminoacid_seq = list(seq.replace('*', '')[aacid_length:])
        # aminoacid_seq.append(is_targeted)
        prepared_seq.append(aminoacid_seq)
    return prepared_seq


def import_dicty_seq_from_file(records):
    dicty_db_entry = []
    for rec in records:
        # Extract the gen_id , sequence_id. They are essential for the later prediction,so they have to be checked.
        ids = rec.description.split("|", 2)
        if ids[0].strip().startswith("DDB") and ids[1].strip().startswith("DDB"):
            sequence_id = ids[0].strip()
            gen_id = ids[1].strip()
        else:
            print("Could not set Entry. Some Error occured for sequence_id and gen_id", rec.description)
            sequence_id = ""
            gen_id = ""
        db_entry = {
            "id": gen_id,
            "gen_id": gen_id,
            "sequence_id": sequence_id,
            "description": rec.description,
            "amino_acid_sequence": str(rec.seq),
            "created_at": datetime.datetime.utcnow(),
            "datasource": "https://testdb.dictybase.org/downloads"
        }
        dicty_db_entry.append(db_entry)
    return dicty_db_entry


def read_in_aminoacids(rel_path, aacid_length, is_targeted):
    raw = list(SeqIO.parse(rel_path, "fasta"))
    return prepare_sequence_for_machine_learning(raw, aacid_length, is_targeted)


def list_to_dataframe(db_proteins):
    proteins = []
    sequences = []
    for db_protein in db_proteins:
        protein = {
            "gen_id": db_protein['gen_id'],
            "sequence_id": db_protein['sequence_id'],
            "description": db_protein['description'],
            "amino_acid_sequence": str(db_protein['amino_acid_sequence']).replace('*', ''),
        }
        proteins.append(protein)
        sequences.append(SeqRecord(Seq(str(db_protein['amino_acid_sequence']).replace('*', ''))))

    df = pd.DataFrame(proteins)
    return df, sequences


def merge_and_convert(neg_targeted, negatives, pos_targeted, positives, min_length_aminoacids):
    seq_concat_data = []
    seq_concat_data.extend(positives)
    seq_concat_data.extend(negatives)
    targeted_concat = []
    targeted_concat.extend(pos_targeted)
    targeted_concat.extend(neg_targeted)

    # Erstellt den Header für die Stellen vom C-terminus (-12, -11, -10, ... -1)
    columns = [str(x) for x in list(range(min_length_aminoacids, 0))]

    # Dataframes mit Einordnung der Aminosäuresequenzen zum Einbuchstabencode
    df_seq = pd.DataFrame(np.array(seq_concat_data), columns=columns)
    df_cat = df_seq.astype(CategoricalDtype(
        categories=aminoacid_onelettercode))
    # Übersetzung der Aminosäuren in binäre Zahlwerte
    seq_dummies = pd.get_dummies(df_cat)
    # Klassifikationen 0/1 der zusammengeführten Sequenzen:
    df_targeted = pd.DataFrame(np.array(targeted_concat), columns=['Targeted'])
    return columns, df_targeted, seq_dummies


def create_svm_with_params(svm, cterminus_length, path_positive_learning_set, path_negative_learning_set):
    negatives, neg_targeted = read_in_aminoacids(path_negative_learning_set,
                                                 cterminus_length, not_peroxisomal)
    positives, pos_targeted = read_in_aminoacids(path_positive_learning_set,
                                                 cterminus_length, peroxisomal)

    # Zusammenführung der peroxisomaler/nicht pox Datensatz.
    columns, df_targeted, seq_dummies = merge_and_convert(neg_targeted, negatives, pos_targeted, positives,
                                                          cterminus_length)

    # Vorbereitung der Lerndaten für die Kreuzvalidierung
    X = seq_dummies.to_numpy()
    y = np.ravel(df_targeted)

    # Trainieren der svm mit allen Lerndaten
    model = svm.fit(X, y)
    return model, columns

### START des SKRIPTS ###

# Festlegen der SVM Parameter und Lerndatensätze
path_positive_filtered_learning_set = "./datasets/pts1_filtered_peroxisomal_sequences.fasta"
path_negative_learning_set = "./datasets/raw_non_peroxisomal_sequences.fasta"
optimal_svm = SVC(kernel='rbf', C=2.0, gamma=0.092, random_state=42)
aa_len = -14

# Training der SVM mit Lerndatensätzen. Erstellung des statistischen Models.
model, columns = create_svm_with_params(optimal_svm, aa_len, path_positive_filtered_learning_set,
                                        path_negative_learning_set)

path_dicty_seq = "./datasets/discoideum_polypeptide.fasta"
dicty_proteins = import_dicty_seq_from_file(list(SeqIO.parse(path_dicty_seq, "fasta")))

# Vorbereitung der Aminosäuresequenzen von Dictyostelium für die SVM
df, sequences = list_to_dataframe(dicty_proteins)
prepared_seq = prepare_dicty_sequences(sequences, aa_len)
df_seq_dicty = pd.DataFrame(np.array(prepared_seq), columns=columns)
df_cat = df_seq_dicty.astype(CategoricalDtype(categories=aminoacid_onelettercode))

# Anwendung der trainierten PTS1-SVM auf die Proteine von Dictyostelium!
dicty_predictions = model.predict(pd.get_dummies(df_cat))

# Vorbereitung zum Speichern der Vorhersagen
s1 = pd.Series(dicty_predictions, name='Targeted')
result = pd.concat([df, df_seq_dicty, s1], axis=1)
predicted_proteins = result.loc[result['Targeted'] == 1]
print("Anzahl vorhergesagter peroxisomaler Proteine mit einer PTS1")
print(len(predicted_proteins))
print(predicted_proteins)
print("Vollständige Übersicht findet sich als Excel unter ./results/predicted_d_discoideum_proteines_with_pts1_file.xlsx")

with pd.ExcelWriter('./results/predicted_d_discoideum_proteines_with_pts1_file.xlsx') as writer:
    predicted_proteins.to_excel(writer, sheet_name='pox_proteine_dicty')
