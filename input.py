import numpy as np
from astropy.io import fits
import astropy.constants as c
import astropy.units as u
import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
import seaborn as sns
from glob import glob
import scipy.ndimage

sns.set_style("white", {"xtick.direction": "in", "ytick.direction": "in"})
    
    

class Page:
    def __init__(self, path):
        self.path = path
        
        self.get_fits()
        # self.set_rms(rms, xy_extent)
        
    def get_fits(self):
        # Get header, image array and rms, convert to micro Jy
        fits_file = fits.open(self.path)
        self.head = fits_file[0].header
        self.im = np.squeeze(fits_file[0].data[0])
        
        # Read in header spatial info to create coordinate grid
        nx = self.head['NAXIS1'];           ny = self.head['NAXIS2']
        xpix = self.head['CRPIX1'];         ypix = self.head['CRPIX2']
        xval = self.head['CRVAL1'];         yval = self.head['CRVAL2']
        xdelt = self.head['CDELT1'];        ydelt = self.head['CDELT2']
        
        self.ra_offset  = np.array( ( (np.arange(nx) - xpix + 1) * xdelt) * 3600)
        self.dec_offset = np.array( ( (np.arange(ny) - ypix + 1) * ydelt) * 3600)
    
        self.wav = (c.c / (self.head['CRVAL3']*u.Hz)).to('mm')
        
    def set_rms(self, rms=None, xy_extent=None):
        try:
            readme_path = self.path[:self.path.rfind('/') + 1] + 'README'
            with open(readme_path, 'r') as f:
                readme = f.read()
                
            rms_ind = readme.find('RMS')
            self.rms = np.float(readme[rms_ind:].split()[1])
            unit = readme[rms_ind:].split()[2]
            
            print('Taking the rms from the following:\n\n')
            print(readme[rms_ind-150:rms_ind+25] + '\n\n')
        except IOError: 
            pass

        if rms is not None:
            self.rms = rms
        elif xy_extent is not None:
            rms_array = np.copy(self.im); 
            x_ind = np.argmin(np.abs(np.abs(self.ra_offset) - xy_extent[0]))
            y_ind = np.argmin(np.abs(np.abs(self.dec_offset) - xy_extent[1]))
            rms_array[x_ind:-x_ind, y_ind:-y_ind] = np.nan
            self.rms = np.sqrt(np.nanmean(rms_array**2)) 
        
        try:
            self.snr = np.array([np.nanmin(self.im), np.nanmax(self.im)]) / self.rms
            print('RMS is {}'.format(self.rms))
            print('SNR ranges between {} and {}'.format(*self.snr))
        except AttributeError:
            self.rms = None
    
    def make_page(self, extent=None, levels=None, channel_ind=None):
        self.fig, self.ax = plt.subplots(figsize=(7,9))
        self.ax.grid(False); self.ax.set_aspect('equal')
        self.ax.set_xlabel(r'$\Delta \alpha$ (")')
        self.ax.set_ylabel(r'$\Delta \delta$ (")')
        
        im = self.im if channel_ind is None else self.im[channel_ind]
        
        if extent != None:
            xmin, xmax = -extent, extent;  ymin, ymax = -extent, extent
            self.ax.set_xlim(xmax, xmin);       self.ax.set_ylim(ymin, ymax)
    
        if levels is not None:
            contour_levels = np.array(levels) * np.nanmax(self.im)
            for_cbar = self.ax.contour(self.ra_offset, self.dec_offset,     
                im, levels=contour_levels, colors='k', 
                linewidths=0.75, extend='both')
        else:
            for_cbar = self.ax.imshow(im,
                extent=[self.ra_offset[0], self.ra_offset[-1], self.dec_offset[-1], self.dec_offset[0]])
                
            # Create the colorbar
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="6%", pad=0.32)
        
        # boundaries = [self.im.min()] + list(contour_levels) + [self.im.max()]
        cbar = self.fig.colorbar(for_cbar, ax=self.ax, cax=cax, 
            ticks=[0], spacing='proportional', orientation='vertical')
        cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), ha='right')
        cbar.ax.yaxis.set_tick_params(pad=10, length=5, direction='out')
        
        try:
            beam_ellipse_color = 'k'
            bmin = self.head['bmin'] * 3600.
            bmaj = self.head['bmaj'] * 3600.
            bpa  = self.head['bpa']

            xmin, xmax = self.ax.get_xlim(); ymin, ymax = self.ax.get_ylim()
            el = mpl.patches.Ellipse(xy=[xmax - (xmax - xmin)/10, ymin + (ymax - ymin)/10], 
                width=bmin, height=bmaj, angle=-bpa,
                edgecolor='k', hatch='///', facecolor='none', zorder=10)
            self.ax.add_artist(el)
        except KeyError: 
            pass
            
    def add_text(self, text_dict):
        for align, string in text_dict.items():
            if align == 'left':
                # x = 1 - page.ax.get_position().extents[2]
                x,y = page.ax.get_position().extents[:2]
            elif align == 'right':
                x,y = page.ax.get_position().extents[2:0:-1]
            self.fig.text(x, y, s=string, horizontalalignment=align)
        
# page = Page('fits_files/AU_Mic_cdaley.fits')
# page.make_page(extent=5, levels=[-0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
# page.add_text({'left' : 'Cail Daley', 'right' : r'AU Mic ALMA {:.2g}'.format(page.wav)})
# plt.savefig('AU_Mic', dpi=100)
# plt.show()
# 
# page = Page('fits_files/49_Ceti_cdaley.fits')
# page.make_page(extent=5, levels=[-0.21, 0.21, 0.4, 0.6, 0.8])
# page.add_text({'left' : 'Jesse Lieman-Sifry & Cail Daley', 'right' : r'49 Ceti ALMA {:.2g}'.format(page.wav)})
# plt.show()
# 
# page = Page('fits_files/TYC_4496_fenclada.fits')
# page.make_page(levels=[-0.1, 0.1, 0.3, 0.5, 0.7, 0.9], extent=3)
# page.add_text({'left' : 'Francisco Encalada', 'right' : 'SMA TYC4496-780-1 {:.2g}'.format(page.wav)})
# plt.show()

# sigma = np.abs(0.33/ (page.head['CDELT1'] * 3600) ) / 2.355
# page=Page('ALMA_data/Elia_2-27_Perez/CLEANcontinuum.sourceElias227.image.fits'); 
# page.im -= 0.87 * scipy.ndimage.gaussian_filter(page.im, sigma=sigma); 
# page.make_page(extent=2, levels=[-0.0192, 0.0192, 0.1, 0.25, 0.5, 0.7, 0.9]); plt.show()

# files = glob('ALMA_data/TW_Hya_Cleeves/*TW_Hya*cube*.pbcor.fits'); files.sort()
# files[2.]
# for file in files:
#     try:
#         page = Page(file)
#     except:
#         pass
#     print(file, page.head['CRVAL3'], page.head['CDELT3'])

# page = Page('ALMA_data/TW_Hya_Cleeves/uid___A001_X88f_X283.TW_Hya_sci.spw29.cube.I.pbcor.fits')
# spectrum = np.nansum(page.im, axis=(1,2)); center_ind = spectrum.argmax() + 1
# page.make_page(levels=[-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], extent=4, channel_ind=center_ind-15); plt.show()
# for i in range(center_ind-15, center_ind+15, 1):
#     page.make_page(levels=[-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], extent=4, channel_ind=i)
#     plt.show()

# #     f
