"""
Plot GNN training metrics across all microscopes
Converted from Jupyter notebook to Python script
"""

# Set matplotlib to use non-interactive backend
import matplotlib
matplotlib.use('Agg')

from sklearn.neural_network import MLPClassifier
import pickle
import os
import pandas as pd
import anndata
import scanpy as sc
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import torch
import math
from sklearn.model_selection import train_test_split, KFold
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp
from sklearn import svm
from sklearn.linear_model import LogisticRegression

from statannotations.Annotator import Annotator
import torch
import pandas as pd
import torch.nn.functional as F
import torch.nn as nn
from torch_geometric.nn import global_mean_pool, global_max_pool, SAGPooling, TopKPooling, global_add_pool
from torch_geometric.nn import DenseGCNConv
import torch_geometric.transforms as T
from torch_geometric.nn import GATConv, GCNConv
from torch_geometric.data import Dataset, DataLoader, DenseDataLoader
import os
import pickle
import random
import math
import numpy as np
from torchmetrics.classification import MulticlassAUROC, MulticlassAccuracy
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import sys
import umap
import matplotlib.pyplot as plt
import seaborn as sns
import time
from datetime import datetime
from joblib import Parallel, delayed

import os.path as osp
from torch_geometric.explain import Explainer, GNNExplainer
# from dgl.nn import AvgPooling, GNNExplainer
from pathlib import Path

sns.set_style("whitegrid")
sns.set(font_scale=2)

# user directories and inputs
# excel file of model paths to compare
dfMicroscopes = pd.read_excel(Path(r"..\data\01 Human lung tissues compared\01 Human lung microscopes compared.xlsx"))
dfMicroscopes.loc[dfMicroscopes['GraphModelPath'].drop_duplicates().index, :] # widefield, localized, and confocal

# Fix paths: if they start with 'data/' instead of '../data/', add the '../' prefix
if not dfMicroscopes['GraphModelPath'].iloc[0].startswith('..'):
    dfMicroscopes['GraphModelPath'] = '../' + dfMicroscopes['GraphModelPath']

print(dfMicroscopes)

# folder to save plot screenshots
screenshotSavePath = Path(r"..\figures")
assert screenshotSavePath.exists()

# Read metrics for all trained GNN models. Select best based on validation ACC, AUC, accuracy, etc.
# read metrics in parallel
def readMetrics(path, microscope):

    dfMetrics = pd.read_csv(path)

    # add model structure info
    structure = path.parts[-6:]

    # select best performing validation score
    valScores = [c for c in dfMetrics.columns if 'val' in c]
    # dfMetrics.dropna(subset = valScores, inplace = True)
    score = dfMetrics[valScores].max() # may need to change this
    score['Structure'] = structure
    # add microscope info
    score['Microscope'] = microscope

    assert not score.isna().any()

    return score

dfAll = []
for row1 in dfMicroscopes.itertuples(): # each microscope type

    # search for all metrics files
    metricFiles = [f for f in Path(row1.GraphModelPath).rglob('*.csv') if 'metrics' in f.stem]

    # read metrics in parallel
    results = Parallel(n_jobs = 10, prefer = 'threads', verbose = 5)\
        (delayed(readMetrics)(path = f, microscope = row1.Microscope) for f in metricFiles)
    dfAll.extend(results)

# concatenate results
dfAll = pd.concat(dfAll, axis=1).T
print(dfAll)

# rename columns for consistency
dfAll.rename(columns = {'val_acc': 'Accuracy',
                        'val_f1': 'F1',
                        'val_auc': 'AUC'}, inplace = True)
# relabel microscopes for consistency
dfAll['Microscope'].replace({'Widefield': 'W',
                             'Localized': 'SRW',
                             'Confocal': 'Airyscan'},
                             inplace = True)

# Plot comapre microscope GNN model performances
def saveFigLabelTime(fig): # save figure with current time as label
    now = datetime.now()
    now = now.strftime('%d%b%Y_%H%M%S')
    fileOut = now + '.png'
    fileOut = os.path.join(screenshotSavePath, fileOut)
    plt.gcf()
    plt.savefig(fileOut, dpi = 300, bbox_inches = 'tight', pad_inches = 0)

    time.sleep(1) # wait to prevent overwrite

    return None

#  reformat metrics
metrics = ['Accuracy', 'F1', 'AUC']
dfSub = dfAll.melt(id_vars = ['Structure', 'Microscope'], value_vars = metrics, value_name='Score', var_name='Metric')
dfSub['Score'] = dfSub['Score'].astype(float)

x = 'Metric'
y = 'Score'
hue = 'Microscope'

fig, ax = plt.subplots(dpi = 300)
sns.barplot(data = dfSub, x = x, y = y, hue = hue, ax = ax)
ax.set_title('GNN Classifier Performance')
ax.legend(bbox_to_anchor=(1, 1), loc='upper left')

# add statistical annotation
scopes = dfSub[hue].unique()
boxPairs = [((m, scopes[0]), (m, scopes[1])) for m in metrics] + [((m, scopes[0]), (m, scopes[-1])) for m in metrics] + \
    [((m, scopes[1]), (m, scopes[-1])) for m in metrics]

annotator = Annotator(ax, boxPairs, data=dfSub, x=x, y=y, hue=hue)
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', comparisons_correction='bonferroni')
annotator.apply_and_annotate()

# save figure
saveFigLabelTime(fig)
plt.close(fig)

print("Script completed successfully!")
