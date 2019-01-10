import methods

# By pass warnings
#====================================================================
import warnings
def warn(*args, **kwargs): pass
warnings.warn = warn

# Define Important Library
#====================================================================
import pandas as pd
import numpy as np

# Start from 1 always, no random state
#====================================================================
np.random.seed(seed=1)

# scikit-learn library import
#====================================================================
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from mlxtend.classifier import EnsembleVoteClassifier

# Necessary Library for Calculation
#====================================================================
from sklearn.metrics import accuracy_score, confusion_matrix, \
    roc_auc_score, average_precision_score, f1_score
#====================================================================
input_file_name = "BenchmarkData.txt"
features_file_name = "Features.csv"

# Here we can generate our features set using different length
#====================================================================
if(1):
    import feature_extractor
    feature_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    feature_extractor.generator(input_file_name, features_file_name, feature_list)

# Load the Featureset:
#====================================================================
D = pd.read_csv(features_file_name, header=None)

# Divide features (X) and classes (y) :
#====================================================================
X = D.iloc[:, :-1].values
y = D.iloc[:, -1].values

# Preprocessing of a large number of data
#====================================================================
from sklearn.preprocessing import Imputer
X[:, 0:X.shape[1]] = Imputer(strategy='mean').fit_transform(X[:, 0:X.shape[1]])

# Encoding y :
#====================================================================
from sklearn.preprocessing import LabelEncoder
y = LabelEncoder().fit_transform(y)



# Define classifiers within a list
#====================================================================
import method_ensemble


# Spliting with 10-FCV :
#====================================================================
from sklearn.model_selection import StratifiedKFold
cv = StratifiedKFold(n_splits=10, shuffle=True)

#====================================================================
from sklearn.pipeline import Pipeline
a = []
b = []
c = []

from methods import add_index
# For SVM
# Feature set 1,2,3,4,5,6,7,8,9,10,11
a = add_index(a,0,16202)


# For LDA
# Feature Set 9, 10, 11, 12, 13, 14, 15, 16, 17
b = add_index(b,0,5)
b = add_index(b,5,8)
b = add_index(b,5546,16384)


# For LR
# Feature Set 12,13,14,15,16,17
c = add_index(c,0,5)
c = add_index(c,16201,16384)

# Making pipeline for all three classifiers and send them into EnsembleVoteClassifier
pipe1 = Pipeline([
               ('sel', method_ensemble.ColumnSelector(a)),
               ('clf', SVC(kernel='rbf', C=4, probability=True, decision_function_shape='ovo', tol=0.1, cache_size=200))])
pipe2 = Pipeline([
               ('sel', method_ensemble.ColumnSelector(b)),
               ('clf', LinearDiscriminantAnalysis(n_components=500))])
pipe3 = Pipeline([
               ('sel', method_ensemble.ColumnSelector(c)),
               ('clf', LogisticRegression(random_state=0, n_jobs=1000))])

#====================================================================
Classifiers = [
    EnsembleVoteClassifier(clfs=[pipe1,pipe2,pipe3], weights=[0.42,0.20,0.38], voting='soft'),
    SVC(kernel='rbf', C=4, probability=True, decision_function_shape='ovo', tol=0.1, cache_size=200),
    LinearDiscriminantAnalysis(n_components=500),
    LogisticRegression(random_state=0, n_jobs=1000),
]
ClassifiersName = [
    'EVC',
    'SVM',
    'LDA',
    'LR',
]

print ('-> Start classification   ...')
# Pick all classifier within the Classifier list and test one by one
#====================================================================
for classifier, cls_name in zip(Classifiers,ClassifiersName):
    np.random.seed(seed=123)
    CM = np.zeros((2, 2), dtype=int)
    CM1 = np.zeros((2, 2), dtype=int)
    fold = 1
    accuracy = []
    auroc = []
    aupr = []
    F1 = []
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    print('___________________________________________________')
    print('Classifier: ' + classifier.__class__.__name__)
    model = classifier

    counter = 1
    for train_index, test_index in cv.split(X, y):
        X_train = X[train_index]
        X_test = X[test_index]

        y_train = y[train_index]
        y_test = y[test_index]

        # Scaling the feature
        from sklearn.preprocessing import StandardScaler, MinMaxScaler
        scale = StandardScaler()
        X_train = scale.fit_transform(X_train)
        X_test = scale.transform(X_test)

        print('Fold: {} '.format(counter))
        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)[:, 1]


        counter += 1
        accuracy.append(accuracy_score(y_pred=y_pred, y_true=y_test))
        auroc.append(roc_auc_score(y_true=y_test, y_score=y_proba))
        aupr.append(average_precision_score(y_true=y_test, y_score=y_proba))
        F1.append(f1_score(y_pred=y_pred, y_true=y_test))

        CM = confusion_matrix(y_pred=y_pred, y_true=y_test)
        CM1 += CM

        TN, FP, FN, TP = CM.ravel()
        print('--------------------------------------------')
        print('| Acc |auROC|auPR | Sp  | Sn  | MCC | F1  |')
        print('--------------------------------------------')
        print('|{:.2f}'.format((accuracy_score(y_pred=y_pred, y_true=y_test) * 100)), end="")
        print('|{:.3f}'.format((roc_auc_score(y_true=y_test, y_score=y_proba))), end="")
        print('|{:.3f}'.format((average_precision_score(y_true=y_test, y_score=y_proba))), end="")

        print('|{:.2f}'.format(((TN / (TN + FP)) * 100)), end="")
        print('|{:.2f}'.format(((TP / (TP + FN)) * 100)), end="")
        print('|{:.3f}'.format(((TP * TN - FP * FN) / (np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))))),
              end="")
        print('|{:.3f}|'.format(np.mean(F1)))
        print('--------------------------------------------')




    print('Average Results:')
    TN, FP, FN, TP = CM1.ravel()
    print('--------------------------------------------')
    print('| Acc |auROC|auPR | Sp  | Sn  | MCC | F1  |')
    print('--------------------------------------------')
    print('|{:.2f}'.format((np.mean(accuracy) * 100)), end="")
    print('|{:.3f}'.format((np.mean(auroc))), end="")
    print('|{:.3f}'.format((np.mean(aupr))), end="")
    print('|{:.2f}'.format(((TN / (TN + FP)) * 100)), end="")
    print('|{:.2f}'.format(((TP / (TP + FN)) * 100)), end="")
    print('|{:.3f}'.format(((TP * TN - FP * FN) / (np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))))), end="")
    print('|{:.3f}|'.format(np.mean(F1)))
    print('--------------------------------------------')
