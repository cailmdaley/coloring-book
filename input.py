%reload_ext autoreload
%autoreload 2
import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
import seaborn as sns
from glob import glob

# try:
#     %matplotlib inline
#     %config InlineBackend.figure_format = 'jpg'
# except:
#     pass
sns.set_style("ticks", {"xtick.direction": "in", "ytick.direction": "in"})
sns.set_context("talk")
    

class Page:
    def __init__(self, path, rms=None, xy_extent=None):
        self.path = path
        self.figsize = (7,7)
        self.get_fits()
        self.set_rms(rms, xy_extent)
        
    def get_fits(self):
        # Get header, image array and rms, convert to micro Jy
        fits_file = fits.open(self.path)
        self.head = fits_file[0].header
        self.im = fits_file[0].data[0][0]
        if np.log10(self.im.max()) >= -3:
            self.unit = 'mJy / beam'
            self.im *= 1e3
        elif np.log10(self.im.max()) >= -6:
            self.unit = r'$\mu$Jy / beam'
            self.im *= 1e6
            
        # Read in header spatial info to create coordinate grid
        nx = self.head['NAXIS1'];           ny = self.head['NAXIS2']
        xpix = self.head['CRPIX1'];         ypix = self.head['CRPIX2']
        xval = self.head['CRVAL1'];         yval = self.head['CRVAL2']
        xdelt = self.head['CDELT1'];        ydelt = self.head['CDELT2']
        
        self.ra_offset  = np.array( ( (np.arange(nx) - xpix + 1) * xdelt) * 3600)
        self.dec_offset = np.array( ( (np.arange(ny) - ypix + 1) * ydelt) * 3600)
    
    def set_rms(self, rms, xy_extent):
        if rms is not None:
            self.rms = rms
        elif xy_extent is not None:
            rms_pixels = []
            for i in range(self.im.shape[0]):
                for j in range(self.im.shape[1]):
                    if np.abs(self.dec_offset[i]) > xy_extent[1] and np.abs(self.ra_offset[i]) > xy_extent[0]:
                        rms_pixels.append(self.im[i,j])
            self.rms = np.sqrt(np.mean(np.array(rms_pixels)**2)) 
    
    def make_page(self, extent=None, sigma_interval=2):
        self.fig, self.ax = plt.subplots(figsize=self.figsize)
        if extent != None:
            xmin, xmax = -extent, extent;  ymin, ymax = -extent, extent
            self.ax.set_xlim(xmax, xmin);       self.ax.set_ylim(ymin, ymax)
            self.ax.grid(False)
    
        # Set x and y major and minor tics
        majorLocator = mpl.ticker.AutoLocator()
        self.ax.xaxis.set_major_locator(majorLocator)
        self.ax.yaxis.set_major_locator(majorLocator)

        # Set x and y labels
        self.ax.set_xlabel(r'$\Delta \alpha$ (")', fontsize=18)
        self.ax.set_ylabel(r'$\Delta \delta$ (")', fontsize=18)
        self.ax.tick_params(which='both', right='on', labelsize=18, direction='in')
        self.ax.tick_params(axis='y', labelright='off', right='on')
            
        try:
            contour_levels = np.arange(sigma_interval, 100, sigma_interval) * self.rms # 3\sigma
            contours = self.ax.contour(self.ra_offset, self.dec_offset, self.im,
                colors='k', levels=contour_levels, 
                linewidths=0.75, linestyles='solid')
            self.ax.contour(self.ra_offset, self.dec_offset, self.im,
                levels=np.flip(contour_levels, axis=0) * -1,
                colors='k', linewidths=0.75, linestyles='dashed')
        except AttributeError:
            cmap = self.ax.imshow(self.im,
                extent=[self.ra_offset[0], self.ra_offset[-1], self.dec_offset[-1], self.dec_offset[0]])

        # Create the colorbar
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(self.ax)
        cax = divider.append_axes("top", size="8%", pad=0.0)
        cmap = mpl.colors.LinearSegmentedColormap.from_list('white', [(1,1,1), (1,1,1)], N=1000)
        cbar = mpl.colorbar.ColorbarBase(cmap=cmap, ax=cax, 
            norm=mpl.colors.Normalize(vmin=self.im.min(), vmax=self.im.max()), 
            orientation='horizontal')
        cbar.ax.xaxis.set_tick_params(which='major', direction='out', length=3,
            bottom='off', top='on', labelsize=8, pad=-2, labeltop='on', 
            labelbottom='off')
        cbar.ax.xaxis.set_tick_params(which='minor', direction='out', length=2, 
            bottom='off', top='on')
        
        try:
            min_sigma = int(self.im.min() / self.rms)
            cbar.set_ticks(list(-contour_levels)[::-1] + [0] + list(contour_levels))
            print([r'${}\sigma$'.format(level) for level in range(-sigma_interval, min_sigma, -sigma_interval)])
            cbar.ax.set_xticklabels(
                [r'${}\sigma$'.format(level) for level in range(-sigma_interval, min_sigma, -sigma_interval)] \
                + ['0'] \
                + [r'${}\sigma$'.format(int(val/self.rms)) for val in contour_levels], 
                fontsize=18)
        except UnboundLocalError:
            pass

        # Overplot the beam ellipse
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
            
    def quickview(self):
            plt.imshow(self.im, origin='lower')
            plt.show(block=False)
            

class ColoringBook:
    # Set seaborn plot styles and color pallete
    def __init__(self):
        self.pages = []
            
    def add_page(self, path):
        self.pages.append(Page(path))
        return self.pages[-1]
        
    def make_pages(self):
        for page in self.pages:
            page.make_page()
            page.quickview()
    
    def add_text(self, subplot_index, x, y, text):
        self.axes[subplot_index].text(x, y, text, fontsize=18,
            path_effects=[mpl.patheffects.withStroke(linewidth=3, 
            foreground="w")])
            
    def add_title(self, title):
        plt.suptitle(self.title)
    
    def show(self):
        plt.show(block=False)
        
    def save(self, path, dpi=400):
        self.fig.savefig(path, dpi=dpi)
        
# files = glob('*.fits')
# files
# f = files[0]

# paths = ['fits_files/' + path for path in ['49_Ceti_cdaley.fits', 'AU_Mic_cdaley.fits', 'TYC_4496_fenclada.fits']]
cb = ColoringBook()
aumic = Page('fits_files/AU_Mic_cdaley.fits', xy_extent=[4,4])
aumic.make_page(extent=4, sigma_interval=5); 
aumic.ax.text(3.6, 3.3, 'AU Mic ALMA 1.4mm')
aumic.ax.text(3.6, 2.7, 'rms = {} '.format(int(aumic.rms.round())) + aumic.unit)
plt.show()

ceti = Page('fits_files/49_Ceti_cdaley.fits', xy_extent=[4,3])
ceti.make_page(sigma_interval=3)
ceti.im.min()/ceti.rms
ceti.ax.text(3.6, 4.1, '49 Ceti ALMA 1.4mm')
ceti.ax.text(3.6, 3.4, 'rms = {} '.format(int(ceti.rms.round())) + ceti.unit)
plt.show()

# tyc = cb.add_page('fits_files/TYC_4496_fenclada.fits')

# ceti.make_page(extent=4.5)
