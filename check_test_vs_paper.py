import anndata
import matplotlib as mp
import pycogaps as pc
import sys

from PyCoGAPS.analysis_functions import *

import logging  as log
log.basicConfig(stream=sys.stderr, level=log.DEBUG)

log.debug('read data')
paper = anndata.read_h5ad('./data/cogapsresult.h5ad')
test = anndata.read_h5ad('./data/testresult_102323.h5ad')

log.debug('create simple plots')
paper_plot = paper.var.plot(kind='line').get_figure()
paper_plot.savefig('paper.png')
test_plot = test.var.plot(kind='line').get_figure()
test_plot.savefig('test.png')

log.debug('create UMAP plots')
plotPatternUMAP(paper, fn='paper')
plotPatternUMAP(test, fn='test')



