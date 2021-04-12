# pts1-prediction-tool
This project offers a classification-algorithm for predicting the peroxisomal targeting signal 1 (PTS1) in a given amino acid sequence (primary structure of a protein).

Current version: 1.0.2
## Installation
``` pip install pts1-prediction-tool```

PyPi: https://pypi.org/project/pts1-prediction-tool/
## Userguide:
Simple example for using the pts1-prediction-tool in your application.
```python
from pts1_prediction_tool.pts1_prediction import PTS1_Predictor

#Instantiates the svm and creates the prediction model
predictor = PTS1_Predictor()

aminoacid_sequence = "MMMMMKLSKMLLLSLSKLSKLSKLSKL"

# Checks a amino acid sequence for an existing PTS1 
result = predictor.check_for_pts1(aminoacid_sequence)

print(result.isPeroxisomal)
```

Run tests with: \
```pipenv shell ``` \
```python -m unittest tests/pts1_prediction_tests.py ```

## Algorithm:
The used classification-algorithm is a support vector machine (svm, sklearn.svm.SVC) from https://scikit-learn.org/. \
This machine learning algorithm was trained to predict the PTS1 in a amino acid sequence (aa_sequence) with a dataset of
514 PTS1/peroxisomal and 11.337 not peroxisomal aa_sequences.
The peroxisomal dataset was generated out of 2324 peroxisomal aa_sequences, which were filtered for the c-terminal-pts1-tripeptide \
(S, A, C, P, H, T, N, Q, E, G, V) / (K, R, H, Q, D, N, S, M) / (L, F, I, M, Y)* 

For the training of the svm, the last 14 C-terminal amino acids of the sequence are used.
The optimal parameters for the svm and the used c-terminal length were determined
by 5-fold-cross validation. For this the perxosiomal and not peroxisomal trainingsets were merged and
separated into 80 % training-sets and 20 % validation-sets\
The final svm has the following statistical average quantities:
* Specifity = 1.0 
* Sensitivity = 0.86
* Precision = 0.98

### Learning data for the svm
The aa_sequences for the learning-sets are downloaded from UniProt (https://www.uniprot.org/), 23.07.2020.
 

An example application will be published on my website, https://olis-lab.de/ \
@Copyright 2021, Oliver Koch 

## Sources:
1. Scikit-learn: Machine Learning in Python, Pedregosa et al., JMLR 12, pp. 2825-2830, 2011.
1. https://scikit-learn.org/, 11.04.2021
1. https://www.uniprot.org/), 23.07.2020
1. https://biopython.org/, 11.04.2021
