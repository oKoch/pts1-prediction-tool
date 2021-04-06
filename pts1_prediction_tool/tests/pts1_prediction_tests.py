import unittest

from pts1_prediction_tool.pts1_prediction import PTS1_Predictor


class Test_pts1_prediction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Test Start")
        cls.predictor = PTS1_Predictor()

    def setUp(self):
        print("Test Start")

    def tearDown(self):
        print("Test Ende")

    def test_predictor_constants(self):
        self.assertEqual(self.predictor.PEROXISOMAL, 1)
        self.assertEqual(self.predictor.NOT_PEROXISOMAL, 0)
        self.assertEqual(len(self.predictor.AMINOACID_ONE_LETTERCODE), 22)

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    unittest.main()

