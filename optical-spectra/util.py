import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
import astropy.units as u
import pandas as pd
from astropy.table import Table

class spectrum1d(object):
    def __init__(self, wave, flux, error=None):
        #Could be np.asarray, but change fit_spectral_lines.py
        self.wave = np.float_(wave)
        self.flux = np.float_(flux)
        self.error = np.float_(error)

def calc_wavelength(header, pixels):
    if ('CTYPE1' in header.keys()) and (len(header['CTYPE1'])>1):
        if header['CTYPE1'] != 'LINEAR': 
            keep_going = input('Attempting to apply a linear solution to non-linear parameters. Continue anyways? (y), n ')
            if keep_going == 'n':
                sys.exit()
    else:
        print('WARNING: No solution type specified, assuming linear')
    CRVAL1 = header['CRVAL1']
    if 'CRPIX1' in header.keys():
        CRPIX1 = header['CRPIX1']
    else:
        CRPIX1 = 0
    if 'CD1_1' in header.keys():
        CD1_1 = header['CD1_1']
    elif 'CDELT1' in header.keys():
        CD1_1 = header['CDELT1']
    elif 'PC1_1' in header.keys():
        CD1_1 = header['PC1_1']
    else:
        sys.exit('Could not identify linear wavelength term (tried CD1_1, CDELT1, PC1_1)')
    
    wavelength = CRVAL1 + CD1_1*(pixels - CRPIX1)
    return wavelength

def make_1d_spec(ifile):
    ofile = fits.open(ifile)
    print(ofile)
    if 'sum_obj' in ifile or '_tel' in ifile or 'MMIRS' in ifile:
        hdr = ofile['FLUX'].header
        print(hdr)
        data = ofile['FLUX'].data
    if 'coadd' in ifile:
        hdr = ofile[1].header
        print(hdr)
        data_tup = ofile[1].data
        data = [],[],[],[],[]
        for i in range(len(data_tup)):
            data[0].append(data_tup[i][0])
            data[1].append(data_tup[i][1])
            data[2].append(data_tup[i][2])
            data[3].append(data_tup[i][3])
            data[4].append(data_tup[i][4])
        print(data[3])
        data = np.array(data)
    if 'sum_obj' not in ifile and '_tel' not in ifile and 'coadd' not in ifile and 'MMIRS' not in ifile:
        hdr = ofile[0].header
        print(hdr)
        data = ofile[0].data
    print(data.shape)
    if 'Bok' in ifile or 'untrim' in ifile or 'test' in ifile or 'GMOS' in ifile:
        #print(data[0,:].shape, data[0,:][0].shape)
        flux = data#[0,:]#[0]
    elif 'obj' in ifile or '_tel' in ifile or 'MMIRS' in ifile and 'Bok' not in ifile:
        flux = data #[0,:]
    elif 'LBT' in ifile:# or 'MMIRS' in ifile:
        flux = data#[0,:]
    elif '_B' in ifile or 'Binospec' in ifile and 'Bok' not in ifile:
        flux = data[0,:]
    elif 'coadd' in ifile:
        flux = data[2]
    else: #if 'untrim' not in ifile and 'Bok' not in ifile and '_B' not in ifile and 'obj' not in ifile and 'coadd' not in ifile and '_tel' not in ifile and 'MMIRS' not in ifile and 'test' not in ifile and 'LBT' not in ifile:
        flux = data[0,:][0]#[1]
        print('the all others case')
    #print(flux.shape)
   # flux = data
    if 'coadd' not in ifile:
        pix = np.arange(len(flux))+1
        #print(flux)
        wl = calc_wavelength(hdr, pix)
    if 'coadd' in ifile:
        wl = data[0]
        print(wl)
    spectrum = spectrum1d(wl, flux)
    return spectrum


def return_header(ifile):
    ofile = fits.open(ifile)
    print(ofile)
    if 'sum_obj' in ifile or '_tel' in ifile or 'MMIRS' in ifile:
        hdr = ofile['FLUX'].header
    if 'coadd' in ifile:
        hdr = ofile[1].header
    if 'sum_obj' not in ifile and '_tel' not in ifile and 'coadd' not in ifile and 'MMIRS' not in ifile:
        hdr = ofile[0].header
    return hdr

def pretty_labels(wave, dist=1100):
    df = pd.DataFrame(wave)
    avg = df.groupby(round(df[0]//dist)).mean()
    return np.array(avg)


def element(name, lines, style, color, row=1, group=False, low=3900, high=6900, dist=1100, ax=plt, start=0.45, spacer=0.65):
    locat = start+(row*spacer)
    ax.vlines(lines, -4, 12, linestyle=style, color=color,  alpha=0.4)
    
    if group==True:
        lines_g = pretty_labels(lines, dist=dist)     
    else: 
        lines_g=lines
    
    for i in range(len(lines_g)):
        if lines_g[i]>=low and lines_g[i]<=high:
            ax.text(lines_g[i], locat, name, ha='center', color=color, rotation=90)
    return

def linelist(ax=plt, spacer=0.45, start=6.05, low=3900, high=6900):
    element_kwargs = {
        "spacer" : spacer,
        "start"  : start,
        "low": low,
        "high": high
    }
    all_elements = [
        element('H', [4102, 4340, 4861, 6563, 10050, 10940., 12820., 15191.2, 15259.9, 15341.1, 15438.2, 15555.7, 15699.9, 15879.8, 16108.6,
                      16406.4, 16805.7, 17361.2, 19444.5, 21660.], '-', 'orchid', row=0, ax=ax, **element_kwargs),  # 3970, 9546, 18173.2, 18751.0,
        element('He', [3889, 4471, 5876, 6678, 7065, 10830], #[5016, 5876, 6678],
                '-', 'darkslateblue', row=0, ax=ax, **element_kwargs),  # 3889, 4471, 4922, , 7281, 10830, 17002, 18685, 19089, 20581
        #element('He II', [4686., 5411], '--', 'C1', row=1, **element_kwargs),
        #element('O II', [3733], ':', 'C3', row=1, **element_kwargs),#[3713, 3727, 3736, 3740, 3749]
        #element('O II', [3630, 3822, 3977, 4235, 4490], ':', 'C1', row=1, **element_kwargs),
        element('[O II]', [3727.1, 3729.8, 7319, 7330], '-.', 'mediumseagreen', row=0, ax=ax, **element_kwargs),
        element('[O III]', [4363, 4959, 5007], '--', 'olivedrab', row=1, ax=ax, **element_kwargs),
        #element('O VI', [3811, 3834], ':', 'C3', row=4, **element_kwargs),
        #element('C II', [6580], '--', 'C2', row=2, **element_kwargs),  # , 723.4
        #element('C III', [4657], '-.', 'C2', row=2, **element_kwargs),  # 4648, 5696, 5826, 6744
        #element('C IV', [5812], ':', 'C2', row=2, **element_kwargs),
        # element('N III', [4640], '-.', 'C5', row=3, ax=ax, **element_kwargs),  # , 8333, 8500.
        # element('N IV', [7123], ':', 'C5', row=3, ax=ax, **element_kwargs),
        # element('N V', [4604, 4620], ':', 'C5', row=4, ax=ax, **element_kwargs),
        #element('[Fe II]', [7155], '-', 'darkslateblue', row=0, ax=ax, **element_kwargs),
        #element('Fe II', [5169, 5002, 5217, 5227, 5260, 7155], ':', 'C9', row=2, **element_kwargs),  # 3683, 3759, 3761, 4348, 5890, 5896
        #element('Fe II', [4924, 5018, 5169], '-', 'C9', row=2, **element_kwargs),
        #element('Ca II', [3934, 3968, 8498, 8542, 8662], '--', 'C4', row=1, **element_kwargs),  # 3706, 3737,
        # 8202, 8249, , 8912, 8927, 9214, 9320, 9321, 9568, 9599, 9891
        #element('Si II', [5972, 6355], '--', 'olive', row=3, **element_kwargs),  # 504.1, 504.9, 505.6, 567.0, 634.7, 635.9, 637.1
        #element('Ti II', [3759, 3761, 3900., 3913, 4300, 4395, 4444, 4468, 4501, 4534, 4550, 4564, 4572, 5966],
        #        ':', 'gray', row=4, **element_kwargs),
        #element('Sc II',
        #        [4306, 4314, 4321, 4374, 4670, 5031, 5240, 5527, 5658, 5684, 6246, 7709, 8098, 8229, 8261,
        #         8371], '--', 'm', row=4, **element_kwargs),
        #element('Sr II', [4078, 4216, 8506, 8689, 8720, 9644, 10040., 10330., 10920.], '--', 'y', row=2, **element_kwargs),
        #element('Ba II', [3892, 4131, 4166, 4525, 4554, 4900, 4934, 4957, 6142, 6497], '--', 'violet',
        #         row=1, group=False, **element_kwargs), #5854, 
        #element('Ba II', [4554, 6142, 6497], '-', 'mediumpurple', row=2, **element_kwargs),
        #element('Mg II', [4391, 4481], '--', 'gold', row=7, **element_kwargs), # , 6347, 6546, 7877, 7896, 8214, 8235, 8735, 8746, 8824, 8835, 9218, 9244, 9328, 9340, 9632
        #element('C',
        #        [4772, 4932, 5052, 5380, 6001, 6006, 6007, 6011, 6013, 6015, 6588, 7113, 7115, 7117, 7861,
        #         8059, 8335, 9061, 9062, 9078, 9089, 9095, 9112, 9406, 9603, 9621, 9658], '--', 'g', row=8, **element_kwargs),
        #element('O', [7774, 8446, 11290], '-', 'darkgreen', row=0, **element_kwargs),  # 7772, 7774, 7775, 8446, 8447, 9263, 9266
        #element('Na', [5890, 5896], '-', '#D81B60', row=0, **element_kwargs),
        # nebular:
        #element('[O I]', [5577, 6300, 6364], ':', 'darkgreen', row=0, dist=1000, ax=ax, **element_kwargs),
        #element('Mg]', [4571], '-.', 'C7', row=0, **element_kwargs),
        #element('[Ca II]', [7293, 7326], ':', 'mediumpurple', row=1, ax=ax, **element_kwargs),
        #element('[Fe II]', [7155], ':', 'C9', row=0, **element_kwargs),
        #element('telluric', [7050], '--', 'gray', row=0, **element_kwargs),
        #element('Mg', [5180], '-', 'C7', row=0, **element_kwargs),
        #element('[Ni II]', [7378], ':', 'olive', row=0, **element_kwargs),
        #element('[C]', [8727], ':', 'C2', row=0, **element_kwargs),
        #element('Fe', [5167, 6359, 7915, 8050, 8070, 8200, 8300, 8340, 8380], '-', 'C9', row=0, group=True, **element_kwargs),
        #element('galaxy', [2798, 3727, 3934, 3969, 4341, 4861, 4959, 5007, 5890, 5896, 6300, 6548, 6563, 6583, 6717, 6731], ':', 'gray', row=0, ax=ax, **element_kwargs),  # 4341, 6300,
        #element('Si', [1203.0], '-', 'gray', row=1, **element_kwargs),
        
        # CORONAL LINES
        element('Fe VII', [3759, 5160, 5722, 6088], '-', 'C9', row=0, group=False, ax=ax, **element_kwargs),
        element('Fe X', [6376], '-', 'C9', row=0, ax=ax, **element_kwargs),
        element('Fe XI', [7894], '-', 'C9', row=0, ax=ax, **element_kwargs),
        element('Fe XIV', [5304], '-', 'C9', row=0, ax=ax, **element_kwargs),
        # element('S XII', )
    ]
