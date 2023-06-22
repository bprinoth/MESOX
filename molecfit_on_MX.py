import numpy as np
import scipy as sc
import tayph.operations as ops
import tayph.tellurics as tel
import pickle
from astropy.io import fits
import copy
from pathlib import Path
import pdb
import matplotlib.pyplot as plt
from tqdm import tqdm
import tayph.util as ut
import pandas as pd
from astropy.time import Time
import glob


def create_s1d_pkl(outpath, file_list, fits_files, response_file):
    
    s1d_headers_blue, s1d_headers_red = [], []

    for file, ff in tqdm.tqdm(zip(file_list, fits_files)):
        
        store = pd.HDFStore(file, 'r+')

        # Read content from hdf data store
        spec_blue   = store['spec_blue']
        spec_red    = store['spec_red']
        store.close()

        # read the response function
        store = pd.HDFStore(response_file, 'r+')
        response_blue = dict(store['response_blue'])
        response_red = dict(store['response_red'])
        store.close()
        

        orders_blue = spec_blue.index.levels[1]
        orders_blue = orders_blue.drop(91)
        orders_red = spec_red.index.levels[1]
        orders_blue = orders_blue[::-1]
        orders_red = orders_red[::-1]

        header = fits.getheader(ff)
        zpt = Time(header["HIERARCH MAROONX TELESCOPE UTCDATE"]+'T00:00:00.000', format='isot', scale='utc')
        obs = Time(header["HIERARCH MAROONX TELESCOPE UTCDATE"]+'T'+header["HIERARCH MAROONX TELESCOPE UTC"], 
            format='isot', scale='utc')

        header["HIERARCH MAROONX TELESCOPE TIME"] = (abs(zpt-obs)).to_value('sec')
        #print(type(header['HIERARCH MAROONX TELESCOPE MJD']))
        header['HIERARCH MAROONX TELESCOPE MJD'] = float(header['HIERARCH MAROONX TELESCOPE MJD'])
        #print(type(header['HIERARCH MAROONX TELESCOPE MJD']))
        
        if 'b_0340' in ff:
            # blue orders
            s1d_headers_blue.append(header)

        elif 'r_0300' in ff:
            # red orders
            s1d_headers_red.append(header)
        else:
            print("I don't know what to do with this... not appending anything.")
            print("something is wrong, I will pdb you out of this...")
            pdb.set_trace()

        # Now we stich the s1d spectra together:

        # first for red
        s1d_flux_stitched_red = []
        s1d_wl_stitched_red = []
        for k in range(len(orders_red)-1):
            o = orders_red[k]
            p = orders_red[k+1]
            flux_larger_order = spec_red['optimal_extraction'][6][p] / response_red[p]
            wl_larger_order = spec_red['wavelengths'][6][p]

            if k==0:
                wl_smaller_order =spec_red['wavelengths'][6][o]
                flux_smaller_order = spec_red['optimal_extraction'][6][o] / response_red[o]
            else: 
                wl_smaller_order = np.asarray(s1d_wl_stitched_red)
                flux_smaller_order = np.asarray(s1d_flux_stitched_red)
            
            # find the overlap reagion in the smaller and larger order:
            wl_overlap_in_smaller_order = wl_smaller_order[wl_smaller_order > np.min(wl_larger_order)]
            wl_overlap_in_larger_order = wl_larger_order[wl_larger_order < np.max(wl_smaller_order)]
            flux_overlap_in_smaller_order = flux_smaller_order[wl_smaller_order > np.min(wl_larger_order)]
            flux_overlap_in_larger_order = flux_larger_order[wl_larger_order <  np.max(wl_smaller_order)]

            # Now interpolate the flux of the larger order in the overlap region onto the wl axis
            # of the overlap region in the smaller order. 
            #pdb.set_trace()
            inp_flux_overlap_region_larger_order = sc.interpolate.interp1d(
                wl_overlap_in_larger_order, flux_overlap_in_larger_order, fill_value='extrapolate', bounds_error=False)(wl_overlap_in_smaller_order)

            mean_flux_overlap_region = 0.5 * (inp_flux_overlap_region_larger_order + flux_overlap_in_smaller_order)
            wl_concatenated = np.concatenate((wl_smaller_order, wl_larger_order[wl_larger_order> np.max(wl_smaller_order)]))
            flux_concatenated = np.concatenate((    flux_smaller_order[wl_smaller_order < np.min(wl_larger_order)],
                                                    mean_flux_overlap_region,
                                                    flux_larger_order[wl_larger_order > np.max(wl_smaller_order)])
            )


            s1d_wl_stitched_red = wl_concatenated
            s1d_flux_stitched_red= flux_concatenated
            
        # and for blue
        s1d_flux_stitched_blue = []
        s1d_wl_stitched_blue = []
        for k in range(len(orders_blue)-1):
            o = orders_blue[k]
            p = orders_blue[k+1]
            flux_larger_order = spec_blue['optimal_extraction'][6][p] / response_blue[p]
            wl_larger_order = spec_blue['wavelengths'][6][p]

            if k==0:
                wl_smaller_order =spec_blue['wavelengths'][6][o]
                flux_smaller_order = spec_blue['optimal_extraction'][6][o] / response_blue[o]
            else: 
                wl_smaller_order = np.asarray(s1d_wl_stitched_blue)
                flux_smaller_order = np.asarray(s1d_flux_stitched_blue)
            
            # find the overlap reagion in the smaller and larger order:
            wl_overlap_in_smaller_order = wl_smaller_order[wl_smaller_order > np.min(wl_larger_order)]
            wl_overlap_in_larger_order = wl_larger_order[wl_larger_order < np.max(wl_smaller_order)]
            flux_overlap_in_smaller_order = flux_smaller_order[wl_smaller_order > np.min(wl_larger_order)]
            flux_overlap_in_larger_order = flux_larger_order[wl_larger_order <  np.max(wl_smaller_order)]

            # Now interpolate the flux of the larger order in the overlap region onto the wl axis
            # of the overlap region in the smaller order. 
            #pdb.set_trace()
            inp_flux_overlap_region_larger_order = sc.interpolate.interp1d(
                wl_overlap_in_larger_order, flux_overlap_in_larger_order, fill_value='extrapolate', bounds_error=False)(wl_overlap_in_smaller_order)

            mean_flux_overlap_region = 0.5*(inp_flux_overlap_region_larger_order + flux_overlap_in_smaller_order)
            wl_concatenated = np.concatenate((wl_smaller_order, wl_larger_order[wl_larger_order> np.max(wl_smaller_order)]))
            flux_concatenated = np.concatenate((    flux_smaller_order[wl_smaller_order < np.min(wl_larger_order)],
                                                    mean_flux_overlap_region,
                                                    flux_larger_order[wl_larger_order > np.max(wl_smaller_order)])
            )


            s1d_wl_stitched_blue = wl_concatenated
            s1d_flux_stitched_blue= flux_concatenated
        

        # now let's write this stuff into pickle file for blue and red

        with open(f'{outpath}_red/s1ds.pkl', 'wb') as f: 
            pickle.dump([s1d_headers_red,s1d_flux_stitched_red,s1d_wl_stitched_red],f)

        f.close()

        with open(f'{outpath}_blue/s1ds.pkl', 'wb') as f: 
            pickle.dump([s1d_headers_blue,s1d_flux_stitched_blue,s1d_wl_stitched_blue],f)

        f.close()


    return 0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# TASK 1: renaming header keywords + creating stitched spectrum



path_to_raw_fits_files = Path('20220403/')
outpath = '/home/bibi/Papers/TRS/data/2022-04-03'
file_list = sorted(glob.glob(path+'*.hd5'))
fits_files = sorted(glob.glob(path+'*.fits'))
response_file = 'MAROON-X_PHOENIX_RESPONSE_CORRECTORDERS.hd5'

create_s1d_pkl(
    outpath=outpath,
    file_list=file_list,
    fits_files=fits_files,
    response_file=response_file
)


# This is it! Your s1d file is created and you can go back to tayph and run molecfit on it.
# Have fun!