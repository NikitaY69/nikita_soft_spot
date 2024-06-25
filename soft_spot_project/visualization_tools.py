import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Circle,Rectangle,FancyArrowPatch,Wedge
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import product, combinations
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
from packaging import version





def render_snapshot(ax, r, d, box, color='gray',edgecolor='black', alpha=1., show_box=True, zorder=0, set_clim=False, climval=[0., 1.], cmap='cividis',cbar=False, kwargs={}):
    patches = []
    ax.axis('off')
    if show_box: ax.axis([0., box[0]+1., 0.-1., box[1]+1.])

    if(isinstance(color, str)):
        for x, y, r in zip(r[:, 0], r[:, 1], d/2.):
            circle = Circle((x, y), r, edgecolor=edgecolor, facecolor=color, linewidth=0.5, alpha=alpha, **kwargs)
            patches.append(circle)
    else:
        for x, y, r in zip(r[:, 0], r[:, 1], d/2.):
            circle = Circle((x, y), r, edgecolor=edgecolor, linewidth=0.5, alpha=alpha, **kwargs)
            patches.append(circle)


    if(isinstance(color, str)):
        p = PatchCollection(patches, match_original=True)
    else:
        p = PatchCollection(patches, cmap=cmap, alpha=alpha)
        p.set_array(color)
        if(set_clim):p.set_clim(climval)


    ax.add_collection(p)
    if show_box:
        if(box.shape[0]==1):
            ax.add_patch(Rectangle((0, 0), width=box,height=box, color='k', fill=False, linewidth=1))
        else:
            ax.add_patch(Rectangle((0, 0), width=box[0],height=box[1], color='k', fill=False, linewidth=1))

    if(cbar):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(p, ax=ax, cax=cax)
        cbar.ax.tick_params(labelsize=7)

    return None



def render_field(ax, r, u, color='k', scale=1., tick_labels=False, labelsize=5, remove_fraction=0, options={}):
    if u.ndim == 1:
        u = np.reshape(u, r.shape)

    if remove_fraction > 0:
        npart = r.shape[0]
        magnitude = u[:, 0]**2 + u[:, 1]**2
        # from small to large
        idx = np.argsort(magnitude)
        idx = np.copy(np.flip(idx))
        idx = idx[:int((1-remove_fraction)*npart)]
        r = r[idx, :]
        u = u[idx, :]

    ax.quiver(r[:, 0], r[:, 1], u[:, 0], u[:, 1], angles='xy',
              scale_units='xy', scale=scale, color=color,headwidth=6, **options)
    if tick_labels:
        L = np.max(r)
        ax.margins(x=0, y=0)
        num_ticks = 4
        locs = np.linspace(0, L, num_ticks)
        ax.set_xticks(locs)
        ax.set_yticks(locs)
        fmt = "%d"
        ax.xaxis.set_major_formatter(FormatStrFormatter(fmt))
        ax.yaxis.set_major_formatter(FormatStrFormatter(fmt))
        ax.tick_params(labelsize=labelsize)
    else:
        ax.margins(x=0, y=0)
        ax.axis('off')

    ax.set_aspect('equal')



