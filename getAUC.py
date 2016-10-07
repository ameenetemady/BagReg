import sys
import csv
import numpy as np
from numpy import genfromtxt
from sklearn import metrics

def getDicRef(strFilename):
    dicRes = {}
    with open(strFilename) as f:
        for line in f:
            dicRes[line.strip()]= True
    return dicRes

strRefFilename=sys.argv[1]
dicRef=getDicRef(strRefFilename)

# read file into my_data
my_data = []
strProtProbsFilename=sys.argv[2]
with open(strProtProbsFilename, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    next(reader, None)
    for row in reader:
        isInRef= row[0] in dicRef
        my_data.append((row[0], float(row[1]) ,int(isInRef)))

# extract columns
y = np.array([ x[2] for x in my_data ], dtype=np.int)
pred = np.array([ x[1] for x in my_data], dtype=np.float)

# calc AUC(roc)
fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=1)
print("AUC(roc):{:f}".format( metrics.auc(fpr, tpr)))

# calc AUC(pr)
precision, recall, thresholds = metrics.precision_recall_curve(y, pred, pos_label=1)
print("AUC(pr):{:f}".format( metrics.auc(recall, precision)))
