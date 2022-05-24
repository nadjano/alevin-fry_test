#!/usr/bin/env python

# import packages
import pyroe 
import sys
import scipy
import pandas as pd

def andata_to_mtx (ad, destination):

    pd.DataFrame(ad.var.index).to_csv(os.path.join(destination, "genes.tsv" ),   sep = "\t", index_col = False)
    pd.DataFrame(ad.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index_col = False)
    # ad.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep = "\t", index_col = True)
    scipy.io.mmwrite(os.path.join(destination, "matrix.mtx"), ad.X)

if __name__ == "__main__":
    #makes andata from USA
    andata = pyroe.load_fry(sys.arg[1])
    #sace andata in matrix market format
    andata_to_mtx(andata, "alevin_fry/mtx/")



    



