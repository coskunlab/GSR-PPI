"""
Plot dot count per cell. Compare untreated vs treated widefield and confocal
Converted from Jupyter notebook to Python script
"""

# Set matplotlib to use non-interactive backend
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import time
import os
import sys
from tqdm import tqdm, trange
from datetime import date
import nd2reader
from joblib import Parallel, delayed
import tifffile as tf
import xml.etree.ElementTree as ET
import napari
import skimage
from skimage.io import imread
from skimage.measure import block_reduce
from skimage.filters import try_all_threshold
import dask.array as da
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
# from dask_ml.preprocessing import MinMaxScaler
import dask
import matplotlib.pyplot as plt
import nd2
from skimage import exposure, restoration
from joblib import Parallel, delayed
import torch
from cellpose import models, core, plot
import cv2
import scipy
import seaborn as sns
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler
from datetime import datetime
from skimage.morphology import closing, square
from skimage.measure import label
import scanpy as sc
import anndata
from anndata import AnnData
import asyncio
from pathlib import Path
from statannotations.Annotator import Annotator
import itertools

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi = 300)

sns.set(font_scale = 2)
sns.set_style('whitegrid')
np.random.seed(0)

# directories and inputs
# folder of individual pair cell dataframes
singlePairPath = Path(r"..\data\HCC827 P7 plate 002\19Jan2024 plate 002 cycle 1 PLA 5 pairs\04 PKL single cell")
assert singlePairPath.exists()
dfSingleInfo = pd.read_excel(r"..\data\HCC827 P7 plate 002\19Jan2024 plate 002 cycle 1 PLA 5 pairs\25Jan2024_wells_info.xlsx")

# folder of combo pair cell dataframes
comboPairPath = Path(r"..\data\HCC827 P7 plate 002\24Jan2024 cycle 1 Keyence\05 PKL single cell")
assert comboPairPath.exists()

# folder of confocal combo cell dataframes
confocalComboPath = Path(r"..\data\HCC827 P7 plate 002\23Jan2024 cycle 1 5 pairs Airyscan re image\02 PKL single cell")
assert confocalComboPath.exists()
dfComboInfo = pd.read_excel(r"..\data\HCC827 P7 plate 002\23Jan2024 cycle 1 5 pairs Airyscan re image\25Jan2024_wells_info.xlsx")

# folder to save plots
screenshotSavePath = Path(r"..\figures")
assert screenshotSavePath.exists()

# Delay execution if necessary
# time.sleep(1 * 60 * 60) # sec

# Reformat and compute single cell widefield
# reformat each FOV dataframe to single cell and concat
def reformatDfSingleCell(path, dfInfo):

    df = pd.read_pickle(path)
    # check if empty
    if df.size == 0:
        return None

    if 'Untreated' in path.stem:
        treatment = 'Untreated'

    elif 'Treated' in path.stem:
        treatment = 'Treated'

    else: # look up treatment from dfSingleInfo
        well = path.stem
        treatment = dfInfo[dfInfo['WellLabel'] == well]
        assert len(treatment) == 1
        treatment = treatment['Drug'].values[0]


    # compute single cell aggregates
    idCols = ['Y', 'X', 'Z', 'FOV', 'MaskCytoLabel', 'MaskNucLabel', 'CellLabel', 'Cycle', 'CellRegion']
    markers = df.drop(columns = idCols).columns
    aggForm = {}
    for ii, markerName in enumerate(markers):
        aggForm[markerName] = 'sum'

    dfCell = df.groupby(['FOV', 'CellLabel']).agg(aggForm).reset_index(drop = False)
    dfCell['Treatment'] = treatment

    return dfCell

# read untreated cells
dfFiles = [f for f in singlePairPath.glob('*.pkl')]
# read in parallel
dfAll = Parallel(n_jobs = 2, prefer = 'threads', verbose = 10)\
(delayed(reformatDfSingleCell)(path = fileName, dfInfo = dfSingleInfo)
 for _, fileName in enumerate(dfFiles))

dfSingle = pd.concat(dfAll)
dfSingle.reset_index(drop = True, inplace = True)
print(dfSingle)

# read treated cells
dfFiles = [f for f in comboPairPath.glob('*.pkl')]
# read in parallel
dfAll = Parallel(n_jobs = 2, prefer = 'threads', verbose = 10)\
(delayed(reformatDfSingleCell)(path = fileName, dfInfo = dfSingleInfo)
 for _, fileName in enumerate(dfFiles))

dfCombo = pd.concat(dfAll)
dfCombo.reset_index(drop = True, inplace = True)
print(dfCombo)

# combine untreated and treated
dfAll = pd.concat([dfSingle, dfCombo])
print(dfAll)

# Boxplot dot counts per cell widefield
def saveFigLabelTime(fig): # save figure with current time as label
    now = datetime.now()
    now = now.strftime('%d%b%Y_%H%M%S')
    fileOut = now + '.png'
    fileOut = Path(screenshotSavePath, fileOut)
    fig.savefig(fileOut, dpi = 300, bbox_inches = 'tight', pad_inches = 0)

    # wait so figures do not overwrite each other
    time.sleep(2)

    return None

# boxplot
dotCols = [col for col in dfAll.columns if '_' in col or 'Combo' in col]
# reformat
widefield = dfAll.melt(id_vars = ['FOV', 'CellLabel', 'Treatment'], value_vars = dotCols, var_name = 'PPI', value_name = 'DotCount')

x = 'PPI'
y = 'DotCount'
hue = 'Treatment'

fig, ax = plt.subplots(dpi = 300)
sns.barplot(data = widefield, x = x, y = y, hue = hue, ax = ax, palette = 'magma')
ax.set_ylabel('Count Expression (AU)')
ax.set_title('Widefield PPI Counts Per Cell')
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend(title = '', bbox_to_anchor = (1, 1), loc = 'upper left')

# add stat test
pairLabels = [f.get_text() for f in ax.get_xticklabels()]
boxPairs = [((p, 'Untreated'), (p, 'Treated')) for p in pairLabels]

annotator = Annotator(ax, boxPairs, data=widefield, x=x, y=y, hue=hue)
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

saveFigLabelTime(fig)
plt.close(fig)

# Plot single cell counts for confocal
# reformat each FOV dataframe to single cell and concat
def reformatDfSingleCell(path, dfInfo):

    df = pd.read_pickle(path)
    # check if empty
    if df.size == 0:
        return None

    # look up treatment from dfSingleInfo
    well = path.stem
    for row1 in dfInfo.itertuples():

        if row1.WellLabel in well:
            treatment = row1.Drug
            break

    # treatment = dfInfo.loc[dfInfo['WellLabel'].isin([well])]
    # print(treatment)
    # assert len(treatment) == 1
    # treatment = treatment['Drug'].values[0]

    # # positive PLA counts
    # df = df.loc[df['Combo'] > 0]

    # compute single cell aggregates
    idCols = ['Y', 'X', 'Z', 'MaskCytoLabel', 'MaskNucLabel', 'Cycle', 'CellRegion']
    # markers = df.drop(columns = idCols).columns
    # aggForm = {}
    # for ii, markerName in enumerate(markers):
    #     aggForm[markerName] = 'sum'
    df.drop(columns = idCols, inplace = True, errors='ignore')
    dfCell = df.groupby(['FOV', 'CellLabel']).sum(numeric_only = True).reset_index(drop = False)

    # dfCell = df.groupby(['FOV', 'CellLabel']).agg(aggForm).reset_index(drop = False)
    dfCell['Treatment'] = treatment

    return dfCell

# read untreated cells
dfFiles = [f for f in confocalComboPath.glob('*.pkl')]
# read in parallel
dfAll = Parallel(n_jobs = 2, prefer = 'threads', verbose = 10)\
(delayed(reformatDfSingleCell)(path = fileName, dfInfo = dfComboInfo)
 for _, fileName in enumerate(dfFiles))

dfConfocal = pd.concat(dfAll)
dfConfocal.reset_index(drop = True, inplace = True)
print(dfConfocal)

# boxplot
dotCols = [col for col in dfConfocal.columns if '_' in col or 'Combo' in col]
# reformat
confocal = dfConfocal.melt(id_vars = ['FOV', 'CellLabel', 'Treatment'], \
                           value_vars = dotCols, var_name = 'PPI', value_name = 'DotCount')

x = 'Treatment'
y = 'DotCount'

fig, ax = plt.subplots(dpi = 300)
sns.barplot(data = confocal, x = x, y = y, ax = ax, palette = 'magma')
ax.set_ylabel('Count Expression (AU)')
ax.set_title('Confocal PPI Counts Per Cell')
# ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
# ax.legend(title = '', bbox_to_anchor = (1, 1), loc = 'upper left')

# add stat test
boxPairs = [tuple(dfConfocal['Treatment'].unique().tolist())]

annotator = Annotator(ax, boxPairs, data=confocal, x=x, y=y)
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

saveFigLabelTime(fig)
plt.close(fig)

# Plot combo confocal and widefield counts
# concat combo confocal and widefield counts
singlePairs = [c for c in widefield.columns if '_' in c]
widefield['Microscope'] = 'Widefield'
confocal['Microscope'] = 'Confocal'
dfCompare = pd.concat([widefield.drop(columns = singlePairs),
                       confocal])
print(dfCompare)

# optional: plot combo pairs only. Drop single pairs
dfCompare = dfCompare.loc[dfCompare['PPI'].str.contains('Combo')]

# boxplot
dotCols = [col for col in dfCompare.columns if '_' in col or 'Combo' in col]

x = 'Microscope'
y = 'DotCount'
hue = 'Treatment'

fig, ax = plt.subplots(dpi = 300)
sns.barplot(data = dfCompare, x = x, y = y, hue = hue, ax = ax, palette = 'magma')
ax.set_ylabel('Count Expression (AU)')
ax.set_title('Multiplex PPI Counts Per Cell')
# ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.legend(title = '', bbox_to_anchor = (1, 1), loc = 'upper left')
ax.set_yscale('log')

# add stat test
pairLabels = [f.get_text() for f in ax.get_xticklabels()]
boxPairs = [((pairLabels[0], treat), (pairLabels[1], treat)) for treat in dfCompare['Treatment'].unique()]

annotator = Annotator(ax, boxPairs, data=dfCompare, x=x, y=y, hue=hue)
annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
annotator.apply_and_annotate()

saveFigLabelTime(fig)
plt.close(fig)

print("Script completed successfully!")
