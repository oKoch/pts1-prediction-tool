import logging

logging.basicConfig(filename='logs.log',level=logging.DEBUG)

class PTS1_Predictor():

    PEROXISOMAL = 1
    NOT_PEROXISOMAL = 0
    AMINOACID_ONE_LETTERCODE = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                                    'T', 'V', 'W', 'Y', 'X', 'U']
    def __init__(self):
        svm = 0


