import anndata
import scanpy as sc
import pandas as pd

pycogapsresult = anndata.read_h5ad("data/pdacresult.h5ad")

# coldata = pd.read_csv("data/PDACcoldata.csv")
# Index(['Unnamed: 0', 'barcode_raw', 'celltype', 'sample_ID',
#        'sample_ID_celltype', 'TN', 'TN_manuscript', 'manuscript', 'nCount_RNA',
#        'nFeature_RNA', 'percent.mt', 'Size_Factor', 'TN_cluster_resolution_5',
#        'TN_assigned_cell_type', 'TN_assigned_cell_type_immune',
#        'TN_assigned_cell_type_immune_specific',
#        'TN_assigned_cell_type_immune_broad', 'cc', 'ccstage',
#        'Classifier_T_duct', 'Classifier_T_Fibroblast_only',
#        'Classifier_T_Fibroblast_Stellate'],
#       dtype='object')
adata = pycogapsresult
# get readable gene names from original object
adata_original = sc.read_h5ad("/Users/jeanette/fertiglab/PDAC_Atlas_Pipeline/PDAC_Peng_Epi.h5ad").T
adata.obs_names = adata_original.obs["gene_short_name"]

# adata = pycogapsresult.T

# adata.obs["cell type"] = list(coldata["TN_assigned_cell_type_immune_specific"])
# adata = pycogapsresult.T
from PyCoGAPS.analysis_functions import *
plotPatternUMAP(adata)

pm = patternMarkers(adata, threshold="cut")

pm["PatternMarkers"]["Pattern1"]
# ['ID3', 'CYP20A1', 'ST3GAL6', 'VGLL4', 'PCDHB14', 'CTSO', 'YARS', 'HMGCLL1', 'ITIH1', 'NCDN', 'TRIM10', 'SMYD2',
# 'HBEGF', 'DHX16', 'TMCO6', 'HSPA1A', 'PHAX', 'EPO', 'CLCN3', 'LSMEM2', 'LY6G5B', 'NDFIP1', 'MICB', 'RUFY1', 'NT5DC2',
# 'ADAMTS4', 'DHX30', 'DVL1', 'PPIL4', 'CBWD2', 'NDST1']

pm["PatternMarkers"]["Pattern2"]
# ['LARP1', 'RIT1', 'SCMH1', 'PCDHA12', 'EGR1', 'COL9A2', 'SMIM13', 'BMP5', 'B3GNT7', 'PLA2G5', 'TPBG', 'EIF3I',
# 'TRIM38', 'MTMR12', 'CAST', 'HMGCS1', 'EYA4']

pm["PatternMarkers"]["Pattern3"]
# ['G3BP1', 'NLRC4', 'ACSL6', 'HSPA4L', 'PRSS56', 'FDCSP', 'FAM8A1', 'GULP1', 'MRPS18B', 'RNF13', 'SNX9', 'CTNND2',
# 'SGCB', 'BDH1', 'AKAP9', 'RWDD4', 'HSD17B4', 'CLIC1', 'BSDC1', 'B4GALT4', 'UBA6', 'NCAPH', 'ECI2', 'MAPK9', 'MLF1',
# 'CEP104', 'ADD1', 'MZT2A', 'LCA5', 'CKAP2L', 'SPRY1', 'UGT8']

pm["PatternMarkers"]["Pattern4"]
# ['EBF1', 'ATP6V1B1', 'CXCL5', 'RWDD3', 'SYF2', 'AP2M1', 'GUCA2B', 'ATP8A1', 'APOBEC2', 'CHDH', 'LRRC1', 'MAP2',
# 'MAST4', 'HMGN3', 'UBE3D', 'TMEM192', 'KLHDC8B', 'TSPYL4', 'IFI44', 'CAPN13', 'C4A', 'EAF1', 'PAK1IP1', 'ORC1',
# 'P2RY12', 'SCN3A', 'TBCCD1', 'GCC2', 'ALDH5A1', 'MBD4', 'FAM162B', 'DROSHA', 'CLDN1', 'ANKRD31', 'NPAS2', 'SNX25',
# 'HADHA', 'SLC25A44', 'PCDHB11']

pm["PatternMarkers"]["Pattern5"]
# ['PDGFRB', 'SPTBN1', 'KCNN2', 'KIF17', 'STAT4', 'KNG1', 'TRANK1', 'VTCN1', 'HNRNPAB', 'CCDC173', 'LSM5', 'CADM2',
# 'BCL10', 'HIST1H4J', 'PRR3', 'FAM193B', 'HARS2', 'PNISR', 'ERI3', 'UBE2E3', 'LTC4S', 'TMEM125', 'GNG12', 'EXTL1',
# 'FAM107A', 'CEP68', 'THRAP3', 'POLE4', 'TRIP6', 'CYP27A1', 'STK19', 'WISP3', 'OSGEPL1', 'ISG20L2', 'ACTG2', 'DHRS3',
# 'SFTPB', 'CLCN2', 'DMP1', 'MGAT5', 'USP19', 'SEPT11', 'AK4', 'ACTL6A', 'RAB1A', 'GALM', 'RARS2', 'SLC25A20', 'YWHAQ',
# 'PTPRF']

pm["PatternMarkers"]["Pattern6"]
# ['LMAN2L', 'ZBTB48', 'ITGAV', 'NRP2', 'MRPS36', 'BLOC1S4', 'KIF4B', 'SLC26A8', 'CYSTM1']

pm["PatternMarkers"]["Pattern7"]
# ['BCL6', 'LAPTM4A', 'TMPRSS11D', 'FAM53A', 'ARSB', 'TRIM46', 'TIMD4', 'FHIT', 'EPHA7', 'BFSP2', 'KHDC1', 'CENPE',
# 'HOOK1', 'BOK', 'STEAP2', 'ITIH4', 'P2RY13', 'GINM1', 'SPP1', 'WLS', 'PPARD', 'CCDC127', 'CYR61', 'RELL2', 'YIPF1',
# 'LOR', 'HAX1', 'LAMTOR2', 'SLC35A4', 'HAT1', 'RGS14', 'GMDS', 'HLTF', 'FAM110D', 'LIPT1', 'RPS27', 'ETV3', 'CSPG5',
# 'SNRNP40', 'LRRC39', 'AMOTL2', 'CBLB', 'AIF1', 'CAD', 'RTKN', 'GCFC2', 'SYNCRIP', 'GTF3C2', 'UBE2J2', 'S100PBP',
# 'HYAL2', 'YBX1', 'DOCK4', 'RABGGTB', 'PARP9']

pm["PatternMarkers"]["Pattern8"]
# ['GFI1', 'KRCC1', 'LRRC8B', 'FGF5', 'HENMT1', 'RNPEP', 'GPBP1', 'PXDC1', 'CLDN11', 'UAP1', 'SAYSD1', 'SPTA1',
# 'IGSF10', 'CREBRF', 'SPINK9', 'ZNF346', 'FAM135A', 'PGM3', 'MANF', 'CENPF', 'PCDHGA9', 'DDX46']