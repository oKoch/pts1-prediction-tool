import logging
import os
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO
from pandas.api.types import CategoricalDtype
from sklearn.svm import SVC

# Sets the right path for the import statements.
CURRENT = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(CURRENT)
sys.path.append(PARENT)
basedir = os.path.abspath(os.path.dirname(__file__))

logging.basicConfig(filename='logs.log', level=logging.DEBUG)

PEROXISOMAL = 1
NOT_PEROXISOMAL = 0
AMINOACID_ONE_LETTERCODE = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                            'T', 'V', 'W', 'Y', 'X', 'U']

# RAW data was taken from UniProt, 23.07.2020
PATH_POSITIVE_FILTERED_LEARNING_SET = os.path.join(basedir, "datasets/pts1_filtered_peroxisomal_sequences.fasta")
PATH_NEGATIVE_LEARNING_SET = os.path.join(basedir, "datasets/raw_non_peroxisomal_sequences.fasta")

'''
This is the length of the C-Terminus that is used for the machine learning model. All aminoacid-sequences are cut to this length
e.g.
N-Terminus-ACDALSDJAPDFJASDJ// CUT //ASJDAPSDAGDHFD-C-Terminus is cut to: ASJDAPSDAGDHFD-C-Terminus
'''
AA_LENGTH = -14


class PTS1_Predictor():

    def __init__(self):
        '''
        The svm-parameters are determined by 5-fold-cross-validation. These offers a svm with the following quality metrics:
        Specifity = 1.0
        Sensitivity = 0.86
        Precision = 0.98
        '''
        svm = SVC(kernel='rbf', C=2.0, gamma=0.092, random_state=42)

        logging.info("Start training the svm machine learning model.")
        self.trained_model, self.columns = self._create_svm_with_params(svm, PATH_POSITIVE_FILTERED_LEARNING_SET,
                                                                        PATH_NEGATIVE_LEARNING_SET)
        logging.info("The svm model is trained and is now usable.")

    def _merge_and_convert(self, neg_targeted, negatives, pos_targeted, positives):
        seq_concat_data = []
        seq_concat_data.extend(positives)
        seq_concat_data.extend(negatives)
        targeted_concat = []
        targeted_concat.extend(pos_targeted)
        targeted_concat.extend(neg_targeted)

        # Erstellt den Header für die Stellen vom C-terminus (-12, -11, -10, ... -1)
        columns = [str(x) for x in list(range(AA_LENGTH, 0))]

        # Dataframes mit Einordnung der Aminosäuresequenzen zum Einbuchstabencode
        df_seq = pd.DataFrame(np.array(seq_concat_data), columns=columns)
        df_cat = df_seq.astype(CategoricalDtype(
            categories=AMINOACID_ONE_LETTERCODE))
        # Übersetzung der Aminosäuren in binäre Zahlwerte
        seq_dummies = pd.get_dummies(df_cat)
        # Klassifikationen 0/1 der zusammengeführten Sequenzen:
        df_targeted = pd.DataFrame(np.array(targeted_concat), columns=['Targeted'])
        return columns, df_targeted, seq_dummies

    def _prepare_sequences_for_machine_learning(self, raw, is_targeted: int):
        aminoacid_sequences = []
        prepared_seq = []
        targeted = []

        for entry in raw:
            # Vereinheitlichung der Aminosäuresequenz zu Großbuchstaben
            seq = str(entry.seq).upper()
            # Aussortieren zu kurzer Aminosäuresequenzen
            if (len(seq) < (-1 * AA_LENGTH)):
                continue

            # Aussortieren der doppelten Aminosäuresequenzen
            if seq not in aminoacid_sequences:
                # Vorbereitung Aussortierung doppelter Sequenzen
                aminoacid_sequences.append(seq)
                # Verkürzen der Sequenz auf letzten stellen vom c-terminus, und umwandlung zu liste, Entfernen des c-terminus sternchen im string
                aminoacid_seq = list(seq.replace('*', '')[AA_LENGTH:])
                targeted.append(is_targeted)
                prepared_seq.append(aminoacid_seq)

        return prepared_seq, targeted

    def _read_in_aminoacid_learning_set(self, rel_path):
        return list(SeqIO.parse(rel_path, "fasta"))

    def _create_svm_with_params(self, svm, path_positive_learning_set, path_negative_learning_set):

        raw_negatives = self._read_in_aminoacid_learning_set(path_negative_learning_set)
        negatives, neg_targeted = self._prepare_sequences_for_machine_learning(raw_negatives, NOT_PEROXISOMAL)

        raw_positives = self._read_in_aminoacid_learning_set(path_positive_learning_set)
        positives, pos_targeted = self._prepare_sequences_for_machine_learning(raw_positives, PEROXISOMAL)

        # Merging and converting of the pos/neg learning datasets
        columns, df_targeted, seq_dummies = self._merge_and_convert(neg_targeted, negatives, pos_targeted, positives)

        X = seq_dummies.to_numpy()
        y = np.ravel(df_targeted)

        # Training of the SVM
        model = svm.fit(X, y)
        return model, columns

    def _prepare_sequences(self, raw_aminoacid_sequence :str):
        prepared_seq = []

        seq = str(raw_aminoacid_sequence).upper()
        if (len(seq) < (-1 * AA_LENGTH)):
            logging.error(f'Aminoacidsequence is to short. The sequence {seq} is not handled.')
            raise Exception(f'Aminoacidsequence is to short. The sequence {seq} is not handled.')

        aminoacid_seq = list(seq.replace('*', '')[AA_LENGTH:])
        prepared_seq.append(aminoacid_seq)

        df_seq_dicty = pd.DataFrame(np.array(prepared_seq), columns=self.columns)
        df_categorical = df_seq_dicty.astype(CategoricalDtype(categories=AMINOACID_ONE_LETTERCODE))

        return pd.get_dummies(df_categorical)

    def check_for_pts1(self, aa_seq: str):
        """
        Checks a aminoacidsequence (aa_seq :str) for an existing PTS1 (peroxisomal targeting sequence 1).
        The aa_seq has to be a fasta like string. The String should start with the N-Terminus and end with the C-Terminus.

        A PTS1PredictionResult is returned. If any error occures, NONE is returned.
        :param aa_seq:
        :return PTS1PredictionResult or None:
        """
        try:
            prepared_aa_sequence = self._prepare_sequences(aa_seq)
            predictions = self.trained_model.predict(prepared_aa_sequence)
            return PTS1PredictionResult(aa_seq=aa_seq, prediction=predictions[0])
        except Exception as ex:
            logging.error(f'The following error occured,{ex}')
            return None


class PTS1PredictionResult():

    def __init__(self, aa_seq: str, prediction: int):
        self.aa_seq = aa_seq
        self.prediction = prediction
        self.isPeroxisomal = False

        if prediction == PEROXISOMAL:
            self.isPeroxisomal = True

# Only for testing purpose:
# if __name__ == "__main__":
#     PTS1_Predictor()