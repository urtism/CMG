import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import datasets
from sklearn.decomposition import PCA
import numpy as np
from pandas import DataFrame, read_csv
import pandas as pd #this is how I usually import pandas
import sys #only needed to determine Python version number
from pyearth import Earth
from sklearn.pipeline import Pipeline
from sklearn.linear_model.logistic import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.feature_extraction import DictVectorizer


def set_missing_values(dataset):
	dataset[dataset=='.']='nan'
	return dataset
	


dataset=r'/home/jarvis/Scrivania/TEST/test_vari/machine_learning/test.dataset'
df = pd.read_csv(dataset,sep='\t')
df= pd.read_table(dataset)

gt_mapping = {
       '0/0': 0,
       '0/1': 1,
       '1/1': 2}

df['GT_GATK'] = df['GT_GATK'].map(gt_mapping)
df['GT_Varscan'] = df['GT_Varscan'].map(gt_mapping)
df['GT_Freebayes'] = df['GT_Freebayes'].map(gt_mapping)

X=df.values[:100,5:]
X=set_missing_values(X)
#print df.columns[12]
y=np.random.randint(2, size=(int(np.shape(X)[0]),))
#print X
#print y 
earth_classifier = Pipeline([('earth', Earth(allow_missing=True)),
                             ('logistic', LogisticRegression())])

#earth_classifier = Pipeline([('earth', Earth(allow_missing=True)),
#                             ('logistic', RandomForestClassifier())])


from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(
	X, y, test_size=0.4, random_state=0)

ec=earth_classifier.fit(X_train, y_train)

y_hat = earth_classifier.predict(X_test)
y_hat_prob = ec.predict_proba(X_test)
acc_hat = ec.score(X_test, y_test)

y_score = ec.fit(X_train, y_train).decision_function(X_test)


i=0
for elem in y_hat:
	print X_test[i,:],y_test[i],elem,y_hat_prob[i],y_score[i]
	i+=1

print acc_hat


import sklearn.metrics as ev
print ev.matthews_corrcoef(y_test,y_hat)
print ev.precision_score(y_test,y_hat)
print ev.recall_score(y_test,y_hat)
print ev.classification_report(y_test,y_hat)
print ev.confusion_matrix(y_test,y_hat)


print ev.roc_curve(y_test,y_hat_prob[:,1])
# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
n_classes = y.shape[0]
#print y_score
for i in range(n_classes):
	fpr[i], tpr[i], _ = ev.roc_curve(y_test, y_hat_prob[:,1])
	roc_auc[i] = ev.auc(fpr[i], tpr[i])



plt.figure()
lw = 2
plt.plot(fpr[2], tpr[2], color='darkorange',
         lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.show()






# Some of these are restricted to the binary classification case:
# matthews_corrcoef(y_true, y_pred[, ...]) 	Compute the Matthews correlation coefficient (MCC) for binary classes
# precision_recall_curve(y_true, probas_pred) 	Compute precision-recall pairs for different probability thresholds
# roc_curve(y_true, y_score[, pos_label, ...]) 	Compute Receiver operating characteristic (ROC)

# Others also work in the multiclass case:
# cohen_kappa_score(y1, y2[, labels, weights]) 	Cohens kappa: a statistic that measures inter-annotator agreement.
# confusion_matrix(y_true, y_pred[, labels, ...]) 	Compute confusion matrix to evaluate the accuracy of a classification
# hinge_loss(y_true, pred_decision[, labels, ...]) 	Average hinge loss (non-regularized)

# Some also work in the multilabel case:
# accuracy_score(y_true, y_pred[, normalize, ...]) 	Accuracy classification score.
# classification_report(y_true, y_pred[, ...]) 	Build a text report showing the main classification metrics
# f1_score(y_true, y_pred[, labels, ...]) 	Compute the F1 score, also known as balanced F-score or F-measure
# fbeta_score(y_true, y_pred, beta[, labels, ...]) 	Compute the F-beta score
# hamming_loss(y_true, y_pred[, labels, ...]) 	Compute the average Hamming loss.
# jaccard_similarity_score(y_true, y_pred[, ...]) 	Jaccard similarity coefficient score
# log_loss(y_true, y_pred[, eps, normalize, ...]) 	Log loss, aka logistic loss or cross-entropy loss.
# precision_recall_fscore_support(y_true, y_pred) 	Compute precision, recall, F-measure and support for each class
# precision_score(y_true, y_pred[, labels, ...]) 	Compute the precision
# recall_score(y_true, y_pred[, labels, ...]) 	Compute the recall
# zero_one_loss(y_true, y_pred[, normalize, ...]) 	Zero-one classification loss.

# And some work with binary and multilabel (but not multiclass) problems:
# average_precision_score(y_true, y_score[, ...]) 	Compute average precision (AP) from prediction scores
# roc_auc_score(y_true, y_score[, average, ...]) 	Compute Area Under the Curve (AUC) from prediction scores
#f = open("filename.txt")
#f.readline()  # skip the header
#data = np.loadtxt(f)

#print data
# import dataset IRIS di priva
iris = datasets.load_iris()
X = iris.data
Y = iris.target

#print iris.DESCR #descrizione del dataset preimpostato
#print iris
#print X,X[:, 0:2]


#CROSS CORRELAZIONE
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(
	X, Y, test_size=0.4, random_state=0)

#print X_train.shape, y_train.shape
#print X_test.shape, y_test.shape
        

#K-FOLD cross validation, train e test hanno una distribuzione delle classi CASUALE
#Each fold is constituted by two arrays: the first one is related to the training set, and the second one to the test set. Thus, one can create the training/test sets using numpy indexing:
# from sklearn.model_selection import KFold

# kf = KFold(n_splits=4)
# for train, test in kf.split(X):
# 	print("%s %s" % (train, test))

# #K-FOLD cross validation stratificata, train e test hanno STESSA distribuzione delle classi
# from sklearn.model_selection import StratifiedKFold

# skf = StratifiedKFold(n_splits=4)
# for train, test in skf.split(X, Y):
# 	print("%s %s" % (train, test))

# X_train, X_test, y_train, y_test = X[train], X[test], Y[train], Y[test]


# x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
# y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5
