"""
Plot PPI counts. Compare HCC827 cells from widefield, localized widefield, and SRRF confocal
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
import dask
import matplotlib
import matplotlib.pyplot as plt
import nd2
from skimage import exposure, restoration
from joblib import Parallel, delayed
import torch
from cellpose import models, core, plot
import cv2
import scipy
from skimage.morphology import closing, square
from skimage.measure import label
import seaborn as sns
import math
from pathlib import Path
import sklearn
from sklearn.neighbors import KDTree
import networkx
from networkx.algorithms.components.connected import connected_components
import h5py
from datetime import datetime

sns.set(font_scale=2)

# user directories and inputs
dfPaths = pd.read_excel(Path(r"../data/HCC827 P6 plate 003/01 HCC827 plate 003 microscopes.xlsx"))

# Fix paths: if they start with 'data/' instead of '../data/', add the '../' prefix
if not dfPaths['PklPath'].iloc[0].startswith('..'):
    dfPaths['PklPath'] = '../' + dfPaths['PklPath']

# well info
dfInfo = pd.read_excel(Path(r"../data/HCC827 P6 plate 003/15Jul2024 cycle 1 Nicky widefield/23Jul2024_wells_info.xlsx"))

screenshotSavePath = Path(r"..\figures")
assert screenshotSavePath.exists()

# Delay execution (optional)
# time.sleep(2 * 60 * 60) # sec

# compute single cell PPI counts
def reformatSingleCell(path, modality, dfInfo):

    df = pd.read_pickle(path)

    # sum dots per cell
    cols = [c for c in df.columns if 'Combo5' in c and 'Intensity' not in c]
    assert len(cols) == 1
    # rename for consistency
    df.rename(columns = {cols[0]: 'Combo5 Dots'}, inplace = True)
    dfCell = df.groupby(['CellLabel']).agg({'Combo5 Dots': 'sum'})
    dfCell.reset_index(drop = False, inplace = True)

    # assign microscope
    dfCell['Microscope'] = modality

    # assign well label
    wellLabel = path.stem
    row = dfInfo.loc[dfInfo['WellLabel'] == wellLabel]
    assert len(row) == 1
    dfCell['WellLabel'] = row['WellLabel'].iloc[0]
    dfCell['Drug'] = row['Drug'].iloc[0]

    return dfCell

dfAll = []
for row1 in dfPaths.itertuples(): # each microscope type

    # get Pkl files
    pklFiles = [f for f in Path(row1.PklPath).glob('*.pkl')]
    # compute single cell in parallel
    results = Parallel(n_jobs=-1, prefer = 'threads', verbose = 10)(delayed(reformatSingleCell)\
                                 (fileName, modality = row1.Microscope, dfInfo = dfInfo) for fileName in pklFiles)

    dfAll.extend(results)

dfCell = pd.concat(dfAll)
# replace microscope labels for consistency
dfCell['Microscope'].replace({'Widefield': 'W',
                              'Localized': 'SRW',
                              'Confocal': 'SRRF'}, inplace = True)
print(dfCell)

# Print number of cells for each microscope group
print(dfCell.groupby(['Microscope']).nunique())

# Plot to compare PPI counts per cell
def saveFigLabelTime(fig): # save figure with current time as label
    now = datetime.now()
    now = now.strftime('%d%b%Y_%H%M%S')
    fileOut = now + '.png'
    fileOut = os.path.join(screenshotSavePath, fileOut)
    fig.savefig(fileOut, dpi = 300, bbox_inches = 'tight', pad_inches = 0)

    # wait so figures do not overwrite each other
    time.sleep(2)
    return None

x = 'Drug'
y = 'Combo5 Dots'
hue = 'Microscope'

fig, ax = plt.subplots(dpi = 300)
ax.set_title('PPI counts per cell')
sns.barplot(x = x, y = y, hue = hue, data = dfCell, ax = ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = 'right')
ax.set_yscale('log')
ax.legend(loc = 'upper left', bbox_to_anchor = (1, 1))

saveFigLabelTime(fig)
plt.close(fig)

# Quantify cell confluency as survival. Mask area / image area
def reformatSingleCell(path, modality, dfInfo, dimY, dimX):

    df = pd.read_pickle(path)

    # compute total mask area / image dimensions
    imageArea = dimY * dimX
    confluency = df['MaskCytoLabel'].size / imageArea # fraction

    dfOut = pd.DataFrame()
    dfOut['Confluency'] = [confluency]

    # assign microscope
    dfOut['Microscope'] = [modality]
    # assign well label
    wellLabel = path.stem
    row = dfInfo.loc[dfInfo['WellLabel'] == wellLabel]
    assert len(row) == 1
    dfOut['WellLabel'] = [row['WellLabel'].iloc[0]]
    dfOut['Drug'] = [row['Drug'].iloc[0]]

    return dfOut

dfAll = []
for row1 in dfPaths.itertuples(): # each microscope type

    # get Pkl files
    pklFiles = [f for f in Path(row1.PklPath).glob('*.pkl')]
    # compute single cell in parallel
    results = Parallel(n_jobs=-1, prefer = 'threads', verbose = 10)(delayed(reformatSingleCell)\
                                 (fileName, modality = row1.Microscope, dfInfo = dfInfo, dimY = row1.DimY, dimX = row1.DimX)
                                 for fileName in pklFiles)

    dfAll.extend(results)

dfFov = pd.concat(dfAll)
# replace microscope labels for consistency
dfFov['Microscope'].replace({'Widefield': 'W',
                              'Localized': 'SRW',
                              'Confocal': 'SRRF'}, inplace = True)
print(dfFov)

# Plot cell confluency (survival) compared across drug groups
x = 'Drug'
y = 'Confluency'
hue = 'Microscope'

fig, ax = plt.subplots(dpi = 300)
ax.set_title('Cell Survival')
ax.set_yscale('log')

sns.barplot(x = x, y = y, hue = hue, data = dfFov, ax = ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, ha = 'right')
ax.legend(loc = 'upper left', bbox_to_anchor = (1, 1))

ax.set_ylabel('Confluency (Survival) Fraction')

saveFigLabelTime(fig)
plt.close(fig)

print("Script completed successfully!")
