import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import cm
import seaborn as sns
sns.set_style("whitegrid")
matplotlib.rc('font',family='DejaVu Sans')
import struct

outdir = "cyl_Re80.out/"

frames = 1001

def load_data(fr):
    fh = "fr_%04d"%fr
    contents = np.loadtxt(outdir+fh+".txt")
    h = contents[0,:]
    nx = int(h[1])
    ny = int(h[2])
    p = contents[1:nx*ny+1,0].reshape(ny,nx)
    u = contents[1:nx*ny+1,1].reshape(ny,nx)
    v = contents[1:nx*ny+1,2].reshape(ny,nx)
    s = np.sqrt(u*u+v*v)
    #vorticity
    w = np.gradient(v,axis=1)+np.gradient(u,axis=0)
    w[0,:]=np.nan
    w[:,0]=np.nan
    w[-1,:]=np.nan
    w[:,-1]=np.nan
    return s,p,w

def get_lims(s,p,w):
    bar = np.zeros_like(p)
    bar[p==0] = 1
    bar[bar==0] = np.nan
    smin = np.min(s)
    smax = np.max(s)
    pmin = np.min(p[p>0])
    pmax = np.max(p)
    wmin = np.nanmin(w)
    wmax = np.nanmax(w)
    return smin,smax,pmin,pmax,wmin,wmax,bar


# first open last frame to get colorbar limits and barrier
s,p,w = load_data(frames-1)
smin,smax,pmin,pmax,wmin,wmax,bar = get_lims(s,p,w)

for fr in range(frames):    
    fh = "fr_%04d"%fr
    s,p,w = load_data(fr)

    # speed plot
    fig = Figure(figsize=(12,16))
    canvas = FigureCanvas(fig)
    ax1 = fig.add_subplot(311)

    cax1 = ax1.imshow(s,cmap=cm.magma,interpolation='bilinear')
    ax1.imshow(bar,cmap=cm.nipy_spectral_r)
    cax1.set_clim(0,np.round(smax,3))
    cbar1 = fig.colorbar(cax1, ticks=[0,np.round(smax/2,3),np.round(smax,3)], orientation='horizontal')
    cbar1.ax.tick_params(labelsize=28)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title("Speed",size=28)

    # density plot
    ax2 = fig.add_subplot(312)
    ax2.contour(p,colors='w',origin='lower',lw=1,alpha=0.5,levels=np.linspace(pmin,pmin+0.045,10))
    cax2 = ax2.imshow(p,cmap=cm.gnuplot,interpolation=None)
    ax2.imshow(bar,cmap=cm.nipy_spectral_r)
    cax2.set_clim(np.round(pmin,3),np.round(pmax,3))
    cbar2 = fig.colorbar(cax2, ticks=[np.round(pmin,3),np.round((pmax-pmin)/2+pmin,3),np.round(pmax,3)], orientation='horizontal')
    cbar2.ax.tick_params(labelsize=28)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title("Density",size=28)

    # vorticity plot
    ax3 = fig.add_subplot(313)

    cax3 = ax3.imshow(w,cmap=cm.plasma,interpolation='bilinear')
    ax3.imshow(bar,cmap=cm.nipy_spectral_r)
    wmax = max([np.abs(wmin),wmax])
    wmin = -wmax
    cax3.set_clim(np.round(wmin/3.,3),np.round(wmax/3.,3))
    cbar3 = fig.colorbar(cax3, ticks=[np.round(wmin/3,3),np.round((wmax/3-wmin/3)/2+wmin/3,3),np.round(wmax/3,3)], orientation='horizontal')
    cbar3.ax.tick_params(labelsize=28)
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_title("Vorticity",size=28)

    fig.subplots_adjust(hspace=0.3)

    fig.tight_layout()
    fig.savefig(outdir+fh+".png",dpi='figure')
    print "Frame",fr
