from astropy.io import fits
from photutils.detection import DAOStarFinder
from photutils.aperture import aperture_photometry, CircularAperture
import numpy as np
import os
from sys import argv
import math

def get_dark(fits_files):
    darks = []
    for file in fits_files:
        if 'dark' in file:
            try:
                with fits.open(file) as hdul:
                    data = hdul[0].data
                    darks.append(data)
            except FileNotFoundError as nf:
                print(nf)
            except Exception as error:
                print(error)
    
    if len(darks) > 0:
        central_dark = darks[0].astype(np.float64)
        for dark in darks[1:]:
            central_dark += dark.astype(np.float64)
        central_dark /= len(darks)
        return central_dark
    else:
        print("No dark frames found.")
        return None

def get_signal(fits_files, dark_frame):
    light_frame = None
    for file in fits_files:
        if 'test_uncal_0' in file:
            light_frame = file
            break
    
    if light_frame is None:
        print("No light frame found.")
        return None
    
    with fits.open(light_frame) as hdul:
        data = hdul[0].data
    
    # calibrate frame
    data = data - dark_frame
    
    # find region of intrest
    std = np.std(data)
    mean = np.mean(data)
    daofind = DAOStarFinder(fwhm=20.0, threshold=5.0 * std)
    stars = daofind(data - mean)
    
    if stars is None:
        print("No stars found.")
        return None

    pos = np.transpose((stars['xcentroid'], stars["ycentroid"]))
    apertures = CircularAperture(pos, r=5.0)
    phot = aperture_photometry(data, apertures)
    
    roi_mean_count = phot['aperture_sum'] / apertures.area
    signal = roi_mean_count.sum() / len(roi_mean_count)
    return signal

def get_readNoise():
    try:
        with fits.open("C:\\Users\\tekka\\Desktop\\OBS\\SNR\\fits-emulator\\rntest.fits") as hdul:
            darkrn = hdul[0].data
        readnoise = np.var(darkrn)
        return readnoise
    except FileNotFoundError as nf:
        print(nf)
        return None
    except Exception as error:
        print(error)
        return None

def calc_SNR(source_signal, dark_signal, read_noise):
    if source_signal is None or dark_signal is None or read_noise is None:
        return None
    snr = source_signal / math.sqrt(source_signal + dark_signal + read_noise)
    return snr

if __name__ == "__main__":
    if len(argv) < 2:
        exit(1)

    folder_path = argv[1]
    if not folder_path.endswith('/'):
        folder_path += '/'
    folder = os.listdir(folder_path)
    fits_files = [file for file in folder if file.endswith('.fits')]
    fits_files = [os.path.join(folder_path, file) for file in fits_files]

    central_dark = get_dark(fits_files)
    if central_dark is not None:
        dark_signal = np.std(central_dark)
    else:
        dark_signal = None
    signal = get_signal(fits_files, central_dark)
    readnoise = 50

    snr = calc_SNR(signal, dark_signal**2, readnoise**2)
    if snr is not None:
        print(f"SNR: {snr}")
    else:
        print("error.")
