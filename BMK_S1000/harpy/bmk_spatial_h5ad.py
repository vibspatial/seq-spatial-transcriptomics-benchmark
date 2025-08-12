import os
import gzip


import numpy as np
import pandas as pd
import cv2
import json
import scanpy as sc


from optparse import OptionParser

parser=OptionParser(description='Generates Spatial h5ad')
parser.add_option('-i', '--input', '--in', default='.', help="BSTViewer_project/subdata/L13_heAuto/ [default: %default]")
parser.add_option('--png', help="BSTViewer_project/he_roi_small.png")
parser.add_option('--library', help="library id")
parser.add_option('-o', '--output', '--out', default='.', help="where to output the gene count matrix [default: %default]")
(opts, args)=parser.parse_args()

print(opts)

if not os.path.exists(opts.output):
        os.mkdir(opts.output)

adata = sc.read_10x_mtx(opts.input, var_names='gene_symbols',cache=True)
barcode_pos = pd.read_csv(f'{opts.input}/barcodes_pos.tsv.gz', compression='gzip', sep='\t', names=["barcode","array_col","array_row"],header=None)
barcode_pos['in_tissue']=1
barcode_pos=barcode_pos[['barcode','in_tissue','array_row','array_col']]
adata.obs = barcode_pos
adata.obs.index = adata.obs['barcode'].to_list()
adata.obs = adata.obs[['in_tissue','array_row','array_col']]

obsm=barcode_pos[['array_col','array_row']]
obsm= obsm.to_numpy()
adata.obsm["spatial"] = obsm

he_img = cv2.imread(opts.png)
he_img=he_img/255
library_id = opts.library
adata.uns["spatial"] = {library_id: {}}
adata.uns["spatial"][library_id]["images"] = {}
adata.uns['spatial'][library_id]['images']['hires'] = he_img.astype(np.float32)
adata.uns['spatial'][library_id]['use_quality']='hires'



def cal_zoom_rate(width, height):
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * np.sqrt(3) / 2.0)
    if std_width / std_height > width / height:
        scale = width / std_width
    else:
        scale = height / std_height
    return scale

zoom_scale = cal_zoom_rate(he_img.shape[1], he_img.shape[0])
adata.uns['spatial'][library_id]['scalefactors']= {"spot_diameter_fullres": zoom_scale, "tissue_hires_scalef": zoom_scale,
                                                   "fiducial_diameter_fullres": zoom_scale, "tissue_lowres_scalef": zoom_scale}
sc.pl.spatial(adata, img_key="hires")

adata.write(f'{opts.output}/adata.h5ad')