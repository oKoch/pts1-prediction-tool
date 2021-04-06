import unittest

from pts1_prediction_tool import pts1_prediction


class Test_pts1_prediction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Test Start")
        cls.predictor = pts1_prediction.PTS1_Predictor()

    def setUp(self):
        print("Test Start")

    def tearDown(self):
        print("Test Ende")

    def test_predictor_constants(self):
        self.assertEqual(pts1_prediction.PEROXISOMAL, 1)
        self.assertEqual(pts1_prediction.NOT_PEROXISOMAL, 0)
        self.assertEqual(len(pts1_prediction.AMINOACID_ONE_LETTERCODE), 22)

    def test_single_prediction(self):
        test_aa_seq = str('MALNFKDKVVIVTGAGGGIGKVYALEFAKRGAKVVVNDLGGSHTGQGSSSKAADKVVEEIKAAGGTAVANYDSVEDGEKIVQTAMDSFGGVDILINNAGILRDVSFG'+
        'KMTDGDWDLVYRVHAKGAYKLSRAAWNHMREKNFGRIIMTSSAAGLYGNFGQANYGSMKMALVGLSNTLAQEGKSKNIHCNTIAPIAASRLTESVMPPEILEQMKPD' +
        'YIVPLVLYLCHQDTTETGGVFEVGAGWVSKVRLQRSAGVYMKDLTPEKIKDNWAQIESFDNPSYPTSASESVSGILAAVNSKPADGESVLVRPPKVAVPKALAATPS' +
        'GSVVVDGYNASKIFTTIQGNIGAKGAELVKKINGIYLINIKKGTNTQAWALDLKNGSGSIVVGAGSTKPNVTITVSDEDFVDIMTGKLNAQSAFTKGKLKISGNMGLA' +
        'TKLGALMQGSKL')

        pred_result = self.predictor.single_prediction(test_aa_seq)
        self.assertEqual(pred_result.prediction, 1)
        self.assertEqual(pred_result.aa_seq, test_aa_seq)


if __name__ == '__main__':
    unittest.main()
