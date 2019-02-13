import numpy as np
import matplotlib.pyplot as plt
from SimPEG import Utils
from scipy.interpolate import (
    LinearNDInterpolator, NearestNDInterpolator
)
import tarfile

# Reproducible science
seed = 12345
np.random.seed(seed)


def download_and_unzip_data(
    url="https://storage.googleapis.com/simpeg/bookpurnong/bookpurnong_inversion.tar.gz"
):
    """
    Download the data from the storage bucket, unzip the tar file, return
    the directory where the data are
    """
    # download the data
    downloads = Utils.download(url)

    # directory where the downloaded files are
    directory = downloads.split(".")[0]

    # unzip the tarfile
    tar = tarfile.open(downloads, "r")
    tar.extractall()
    tar.close()

    return downloads, directory

def plot2Ddata_categorical(
    xyz, data, vec=False, nx=50, ny=50,
    ax=None, mask=None, level=False, figname=None,
    ncontour=10, dataloc=False, contourOpts={},
    levelOpts={}, scale="linear", clim=None
):
    """

        Take unstructured xy points, interpolate, then plot in 2D

        :param numpy.array xyz: data locations
        :param numpy.array data: data values
        :param bool vec: plot streamplot?
        :param float nx: number of x grid locations
        :param float ny: number of y grid locations
        :param matplotlib.axes ax: axes
        :param numpy.array mask: mask for the array
        :param boolean level: boolean to plot (or not)
                                :meth:`matplotlib.pyplot.contour`
        :param string figname: figure name
        :param float ncontour: number of :meth:`matplotlib.pyplot.contourf`
                                contours
        :param bool dataloc: plot the data locations
        :param dict controuOpts: :meth:`matplotlib.pyplot.contourf` options
        :param dict levelOpts: :meth:`matplotlib.pyplot.contour` options
        :param numpy.array clim: colorbar limits

    """
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)

    xmin, xmax = xyz[:, 0].min(), xyz[:, 0].max()
    ymin, ymax = xyz[:, 1].min(), xyz[:, 1].max()
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(x, y)
    xy = np.c_[X.flatten(), Y.flatten()]
    if vec is False:
        F = NearestNDInterpolator(xyz[:, :2], data)
        DATA = F(xy)
        DATA = DATA.reshape(X.shape)
        if scale == "log":
            DATA = np.log10(abs(DATA))

        # Levels definitions
        dataselection = np.logical_and(
            ~np.isnan(DATA),
            np.abs(DATA) != np.inf
            )
        if clim is None:
            vmin = DATA[dataselection].min()
            vmax = DATA[dataselection].max()
        else:
            vmin = np.min(clim)
            vmax = np.max(clim)
            if scale == "log":
                vmin = np.log10(vmin)
                vmax = np.log10(vmax)
        if np.logical_or(
            np.logical_or(
                np.isnan(vmin),
                np.isnan(vmax)
            ),
            np.logical_or(
                np.abs(vmin) == np.inf,
                np.abs(vmax) == np.inf
            )
        ):
            raise Exception(
                """clim must be sctrictly positive in log scale"""
            )
        vstep = np.abs((vmin-vmax)/(ncontour+1))
        levels = np.arange(vmin, vmax+vstep, vstep)
        if DATA[dataselection].min() < levels.min():
                levels = np.r_[DATA[dataselection].min(), levels]
        if DATA[dataselection].max() > levels.max():
                levels = np.r_[levels, DATA[dataselection].max()]

        if mask is not None:
            DATA = np.ma.masked_array(DATA, mask=mask)

        cont = ax.contourf(
            X, Y, DATA, levels=levels,
            vmin=vmin, vmax=vmax,
            **contourOpts
        )
        if level:
            CS = ax.contour(X, Y, DATA, levels=levels, **levelOpts)

    else:
        # Assume size of data is (N,2)
        datax = data[:, 0]
        datay = data[:, 1]
        Fx = NearestNDInterpolator(xyz[:, :2], datax)
        Fy = NearestNDInterpolator(xyz[:, :2], datay)
        DATAx = Fx(xy)
        DATAy = Fy(xy)
        DATA = np.sqrt(DATAx**2+DATAy**2).reshape(X.shape)
        DATAx = DATAx.reshape(X.shape)
        DATAy = DATAy.reshape(X.shape)
        if scale == "log":
            DATA = np.log10(abs(DATA))

        # Levels definitions
        dataselection = np.logical_and(
            ~np.isnan(DATA),
            np.abs(DATA) != np.inf
            )
        if clim is None:
            vmin = DATA[dataselection].min()
            vmax = DATA[dataselection].max()
        else:
            vmin = np.min(clim)
            vmax = np.max(clim)
            if scale == "log":
                vmin = np.log10(vmin)
                vmax = np.log10(vmax)
        if np.logical_or(
            np.logical_or(
                np.isnan(vmin),
                np.isnan(vmax)
            ),
            np.logical_or(
                np.abs(vmin) == np.inf,
                np.abs(vmax) == np.inf
            )
        ):
            raise Exception(
                """clim must be sctrictly positive in log scale"""
            )
        vstep = np.abs((vmin-vmax)/(ncontour+1))
        levels = np.arange(vmin, vmax+vstep, vstep)
        if DATA[dataselection].min() < levels.min():
                levels = np.r_[DATA[dataselection].min(), levels]
        if DATA[dataselection].max() > levels.max():
                levels = np.r_[levels, DATA[dataselection].max()]

        if mask is not None:
            DATA = np.ma.masked_array(DATA, mask=mask)

        cont = ax.contourf(
            X, Y, DATA, levels=levels,
            vmin=vmin, vmax=vmax,
            **contourOpts
        )
        ax.streamplot(X, Y, DATAx, DATAy, color="w")
        if level:
            CS = ax.contour(X, Y, DATA, levels=levels, **levelOpts)

    if dataloc:
        ax.plot(xyz[:, 0], xyz[:, 1], 'k.', ms=2)
    plt.gca().set_aspect('equal', adjustable='box')
    if figname:
        plt.axis("off")
        fig.savefig(figname, dpi=200)
    if level:
        return cont, ax, CS
    else:
        return cont, ax


def plot2Ddata(
    xyz, data, vec=False, nx=50, ny=50,
    ax=None, mask=None, level=False, figname=None,
    ncontour=10, dataloc=False, contourOpts={},
    levelOpts={}, scale="linear", clim=None
):
    """

        Take unstructured xy points, interpolate, then plot in 2D

        :param numpy.array xyz: data locations
        :param numpy.array data: data values
        :param bool vec: plot streamplot?
        :param float nx: number of x grid locations
        :param float ny: number of y grid locations
        :param matplotlib.axes ax: axes
        :param numpy.array mask: mask for the array
        :param boolean level: boolean to plot (or not)
                                :meth:`matplotlib.pyplot.contour`
        :param string figname: figure name
        :param float ncontour: number of :meth:`matplotlib.pyplot.contourf`
                                contours
        :param bool dataloc: plot the data locations
        :param dict controuOpts: :meth:`matplotlib.pyplot.contourf` options
        :param dict levelOpts: :meth:`matplotlib.pyplot.contour` options
        :param numpy.array clim: colorbar limits

    """
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)

    xmin, xmax = xyz[:, 0].min(), xyz[:, 0].max()
    ymin, ymax = xyz[:, 1].min(), xyz[:, 1].max()
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(x, y)
    xy = np.c_[X.flatten(), Y.flatten()]
    if vec is False:
        F = LinearNDInterpolator(xyz[:, :2], data)
        DATA = F(xy)
        DATA = DATA.reshape(X.shape)
        if scale == "log":
            DATA = np.log10(abs(DATA))

        # Levels definitions
        dataselection = np.logical_and(
            ~np.isnan(DATA),
            np.abs(DATA) != np.inf
            )
        if clim is None:
            vmin = DATA[dataselection].min()
            vmax = DATA[dataselection].max()
        else:
            vmin = np.min(clim)
            vmax = np.max(clim)
            if scale == "log":
                vmin = np.log10(vmin)
                vmax = np.log10(vmax)
        if np.logical_or(
            np.logical_or(
                np.isnan(vmin),
                np.isnan(vmax)
            ),
            np.logical_or(
                np.abs(vmin) == np.inf,
                np.abs(vmax) == np.inf
            )
        ):
            raise Exception(
                """clim must be sctrictly positive in log scale"""
            )
        vstep = np.abs((vmin-vmax)/(ncontour+1))
        levels = np.arange(vmin, vmax+vstep, vstep)
        if DATA[dataselection].min() < levels.min():
                levels = np.r_[DATA[dataselection].min(), levels]
        if DATA[dataselection].max() > levels.max():
                levels = np.r_[levels, DATA[dataselection].max()]

        if mask is not None:
            DATA = np.ma.masked_array(DATA, mask=mask)

        cont = ax.contourf(
            X, Y, DATA, levels=levels,
            vmin=vmin, vmax=vmax,
            **contourOpts
        )
        if level:
            CS = ax.contour(X, Y, DATA, levels=levels, **levelOpts)

    else:
        # Assume size of data is (N,2)
        datax = data[:, 0]
        datay = data[:, 1]
        Fx = LinearNDInterpolator(xyz[:, :2], datax)
        Fy = LinearNDInterpolator(xyz[:, :2], datay)
        DATAx = Fx(xy)
        DATAy = Fy(xy)
        DATA = np.sqrt(DATAx**2+DATAy**2).reshape(X.shape)
        DATAx = DATAx.reshape(X.shape)
        DATAy = DATAy.reshape(X.shape)
        if scale == "log":
            DATA = np.log10(abs(DATA))

        # Levels definitions
        dataselection = np.logical_and(
            ~np.isnan(DATA),
            np.abs(DATA) != np.inf
            )
        if clim is None:
            vmin = DATA[dataselection].min()
            vmax = DATA[dataselection].max()
        else:
            vmin = np.min(clim)
            vmax = np.max(clim)
            if scale == "log":
                vmin = np.log10(vmin)
                vmax = np.log10(vmax)
        if np.logical_or(
            np.logical_or(
                np.isnan(vmin),
                np.isnan(vmax)
            ),
            np.logical_or(
                np.abs(vmin) == np.inf,
                np.abs(vmax) == np.inf
            )
        ):
            raise Exception(
                """clim must be sctrictly positive in log scale"""
            )
        vstep = np.abs((vmin-vmax)/(ncontour+1))
        levels = np.arange(vmin, vmax+vstep, vstep)
        if DATA[dataselection].min() < levels.min():
                levels = np.r_[DATA[dataselection].min(), levels]
        if DATA[dataselection].max() > levels.max():
                levels = np.r_[levels, DATA[dataselection].max()]

        if mask is not None:
            DATA = np.ma.masked_array(DATA, mask=mask)

        cont = ax.contourf(
            X, Y, DATA, levels=levels,
            vmin=vmin, vmax=vmax,
            **contourOpts
        )
        ax.streamplot(X, Y, DATAx, DATAy, color="w")
        if level:
            CS = ax.contour(X, Y, DATA, levels=levels, **levelOpts)

    if dataloc:
        ax.plot(xyz[:, 0], xyz[:, 1], 'k.', ms=2)
    plt.gca().set_aspect('equal', adjustable='box')
    if figname:
        plt.axis("off")
        fig.savefig(figname, dpi=200)
    if level:
        return cont, ax, CS
    else:
        return cont, ax


def order_clusters_GM_mean(GMmodel, outputindex=False):
    from SimPEG.Utils import _compute_precision_cholesky
    '''
    order cluster by increasing mean for Gaussian Mixture scikit object
    '''

    indx = np.argsort(GMmodel.means_, axis=0)
    GMmodel.means_ = GMmodel.means_[indx].reshape(GMmodel.means_.shape)
    GMmodel.weights_ = GMmodel.weights_[indx].reshape(GMmodel.weights_.shape)
    if GMmodel.covariance_type == 'tied':
        pass
    else:
        GMmodel.precisions_ = GMmodel.precisions_[
            indx].reshape(GMmodel.precisions_.shape)
        GMmodel.covariances_ = GMmodel.covariances_[
            indx].reshape(GMmodel.covariances_.shape)
    GMmodel.precisions_cholesky_ = _compute_precision_cholesky(
        GMmodel.covariances_, GMmodel.covariance_type)

    if outputindex:
        return indx


def weighted_avg_and_var(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2., weights=weights)
    return (average, variance)
