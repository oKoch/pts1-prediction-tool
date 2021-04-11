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

    def test_prediction_with_peroxisomal_seq(self):
        test_aa_seq = str('MALNFKDKVVIVTGAGGGIGKVYALEFAKRGAKVVVNDLGGSHTGQGSSSKAADKVVEEIKAAGGTAVANYDSVEDGEKIVQTAMDSFGGVDILINNAGILRDVSFG'+
        'KMTDGDWDLVYRVHAKGAYKLSRAAWNHMREKNFGRIIMTSSAAGLYGNFGQANYGSMKMALVGLSNTLAQEGKSKNIHCNTIAPIAASRLTESVMPPEILEQMKPD' +
        'YIVPLVLYLCHQDTTETGGVFEVGAGWVSKVRLQRSAGVYMKDLTPEKIKDNWAQIESFDNPSYPTSASESVSGILAAVNSKPADGESVLVRPPKVAVPKALAATPS' +
        'GSVVVDGYNASKIFTTIQGNIGAKGAELVKKINGIYLINIKKGTNTQAWALDLKNGSGSIVVGAGSTKPNVTITVSDEDFVDIMTGKLNAQSAFTKGKLKISGNMGLA' +
        'TKLGALMQGSKL')

        pred_result = self.predictor.check_for_pts1(test_aa_seq)
        self.assertEqual(pred_result.prediction, 1)
        self.assertEqual(pred_result.aa_seq, test_aa_seq)
        self.assertEqual(pred_result.isPeroxisomal, True)

    def test_prediction_with_not_peroxisomal_seq(self):
        test_aa_seq = str('MALNFKDKVVIVTGAGGGIGKVYALEFAKRGAKVVVNDLGGSHTGQGSSSKAADKVVEEIKAAGGTAVANYDSVEDGEKIVQTAMDSFGGVDILINNAGILRDVSFG'+
        'KMTDGDWDLVYRVHAKGAYKLSRAAWNHMREKNFGRIIMTSSAAGLYGNFGQANYGSMKMALVGLSNTLAQEGKSKNIHCNTIAPIAASRLTESVMPPEILEQMKPD' +
        'YIVPLVLYLCHQDTTETGGVFEVGAGWVSKVRLQRSAGVYMKDLTPEKIKDNWAQIESFDNPSYPTSASESVSGILAAVNSKPADGESVLVRPPKVAVPKALAATPS' +
        'GSVVVDGYNASKIFTTIQGNIGAKGAELVKKINGIYLINIKKGTNTQAWALDLKNGSGSIVVGAGSTKPNVTITVSDEDFVDIMTGMMMMMMMMMMMMMMMMMMMMM')
        pred_result = self.predictor.check_for_pts1(test_aa_seq)

        self.assertEqual(pred_result.prediction, 0)
        self.assertEqual(pred_result.aa_seq, test_aa_seq)
        self.assertEqual(pred_result.isPeroxisomal, False)

    def test_prediction_with_to_short_aa_seq(self):
        test_aa_seq = str('GALMQGSKL')
        pred_result = self.predictor.check_for_pts1(test_aa_seq)
        self.assertIsNone(pred_result)

if __name__ == '__main__':
    unittest.main()
