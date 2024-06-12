'''
To run this file open CASA and define the following variables

dirs = [] # a list of fits file directories to run on
'''

import os
import glob
# target = 'AT2022dbl_'
# dirs = ['S','S1','S2','C','C1','C2','X','X1','X2']

# define an offsourcebox
offsourcebox = '0,0,30,30'

for fits_dir in dirs:
    fits = glob.glob(os.path.join(fits_dir, '*.fits'))
    print(fits)
    for img in fits:
        for i, S in zip([0, 1, 2, 3, 0], ['I', 'Q', 'U', 'V', 'P']):
            fitfile = os.path.join(f'{img.replace(".fits","")}_results.txt')
            estimatesname = 'estimates'

            # get the center x, y coordinates (integer divide the shape from the header)
            npix = 20
            x, y, _, _ = imhead(imagename=img, mode='get', hdkey='shape') // 2
            srcbox = str(x - npix) + ',' + str(y - npix) + ',' + str(x + npix) + ',' + str(y + npix)

            dateobs = imhead(imagename=img, mode='get', hdkey='date-obs')
            bmaj = imhead(imagename=img, mode='get', hdkey='beammajor')['value']
            bmajunit = imhead(imagename=img, mode='get', hdkey='beammajor')['unit']
            bmin = imhead(imagename=img, mode='get', hdkey='beamminor')['value']
            bminunit = imhead(imagename=img, mode='get', hdkey='beamminor')['unit']
            bpa = imhead(imagename=img, mode='get', hdkey='beampa')['value']
            bpaunit = imhead(imagename=img, mode='get', hdkey='beampa')['unit']
            estimates = '0.005,' + str(x) + ',' + str(y) + ',' + str(bmaj) + bmajunit + ',' + str(bmin) + bminunit + ',' + str(bpa) + bpaunit + ',abp'
            estimatesfile = estimatesname + '.txt'
            f = open(estimatesfile, 'w')
            f.write(estimates)
            f.close()

            resultsI = imfit(imagename=img, estimates=estimatesfile, box=srcbox)
            
            resultsI = resultsI['results']['component0']
            xypos = resultsI['pixelcoords']
            q = resultsI['spectrum']['frequency']['m0']['value']
            ra = resultsI['shape']['direction']['m0']['value']
            dra = resultsI['shape']['direction']['error']['longitude']['value']
            dec = resultsI['shape']['direction']['m1']['value']
            ddec = resultsI['shape']['direction']['error']['latitude']['value']

            F = resultsI['flux']['value'][i] * 1000.
            dF = resultsI['flux']['error'][i] * 1000.

            results2 = imstat(imagename=img, box=offsourcebox)['sigma'][0] * 1000.

            # Write results to a file                                                                                                                                                                           
            f = open(fitfile, 'a')
            f.write(dateobs + ', ' + str(q) + ', ' + S + ', ' + str(ra) + ', ' + str(dra) + ', ' + str(dec) + ', ' + str(ddec) + ', ' + str(F) + ', ' + str(dF) + ', ' + str(results2) + '\n')
            f.close()

