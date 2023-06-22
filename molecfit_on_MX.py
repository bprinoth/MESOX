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

def load_order(dp, nr=20, expnr=28):
    """
    dp: datapath to folder with orders and wavelength files
    nr: order number that should be considered
    """
    
    # IMPORT HEADER INFORMATION FROM s1ds.pkl
    with open(dp/'s1ds.pkl', 'rb') as f: 
        s1dheaders, _,_ = pickle.load(f)
        
    
    order = fits.getdata(dp/f'order_{nr}.fits')
    wave = fits.getdata(dp/f'wave_{nr}.fits') * 10. # in Å in vac
    #print(np.shape(order))
    
    order = order[expnr]
    wave = wave[expnr]
    
    p = np.polyfit(wave, order, deg=6)
    poly = np.polyval(p, wave) 
    
  
    return s1dheaders[expnr], order / poly , wave
   
    
def write_file_to_molecfit(molecfit_folder,name,header,wave,spectrum,plot=False, time_average=False):
    spectrum[spectrum<=0]=np.nan
    err = np.sqrt(spectrum)
    
    col1 = fits.Column(name = 'wavelength', format = '1D', array = wave)
    col2 = fits.Column(name = 'flux', format       = '1D', array = spectrum)
    col3 = fits.Column(name = 'err_flux', format   = '1D', array = err)
    
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr = copy.deepcopy(header)
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(molecfit_folder/name,overwrite=True)
    
    return 0

 # add your own paths here!
dp = Path('data/2022-04-03/')
molecfit_input_folder = Path('/data/bibi/Papers/TRS/input_files/2022-04-03')
instrument= 'MAROON-X'
molecfit_prog_folder = Path('/usr/local/molecfit_1.5.9/bin')
python_alias = 'python3'
parname = Path('MAROON-X.par')
parfile = molecfit_input_folder/parname

orders_to_look_at = [27]


for i in tqdm(range(27,28)):

    #print(f'ORDER {i}')
    waves, fxs, transs = [],[],[]
    if i in orders_to_look_at:
       
        # LOAD IN ORDER i, expnr

        header, spectrum, wave = load_order(dp, i)
        
        
        write_file_to_molecfit(molecfit_input_folder,instrument+'.fits',header,wave,spectrum)
        tel.execute_molecfit(molecfit_prog_folder,parfile,molecfit_input_folder,gui=True,
            alias=python_alias)

        for k in range(57):
            header, spectrum, wave = load_order(dp, i, expnr=k)    
            write_file_to_molecfit(molecfit_input_folder,instrument+'.fits',header,wave,spectrum)

            tel.execute_molecfit(molecfit_prog_folder,parfile,molecfit_input_folder,gui=False,
                alias=python_alias)
            wl,fx,trans = tel.retrieve_output_molecfit(molecfit_input_folder/instrument)

            waves.append(wl)
            fxs.append(fx/trans)
            transs.append(trans)
    else:
        for k in range(57):
            header, spectrum, wave = load_order(dp, i, expnr=k)
            wl = wave
            fx = spectrum
            trans = np.ones_like(spectrum)

            waves.append(wl)
            fxs.append(fx/trans)
            transs.append(trans)

        ut.writefits(f'/data/bibi/Papers/TRS/data/2022-04-03/tellurics_wave_{i}.fits', np.asarray(waves))
        ut.writefits(f'/data/bibi/Papers/TRS/data/2022-04-03/tellurics_order_{i}.fits', np.asarray(transs))


# plt.figure()
# for m in range(len(transs)):
#     plt.plot(waves[m], transs[m])

# plt.show()
