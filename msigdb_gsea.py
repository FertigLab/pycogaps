import pickle
import scanpy as sc
import pandas as pd
from PyCoGAPS.analysis_functions import *

with open('./data/spatialhypoxia8patterns.pkl', 'rb') as fp:
    hypoxiaresult8pattern = pickle.load(fp)

hypoxia_hallmarks = pd.read_table("hypoxiagenes", header=None)

genes = list(hypoxia_hallmarks[0])
genes = list(set(hypoxiaresult8pattern.obs_names).intersection(genes))






cogaps_stat = calcCoGAPSStat(hypoxiaresult8pattern, genes)

