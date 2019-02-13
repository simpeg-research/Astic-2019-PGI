from SimPEG.EM.Static import DC, Utils as DCUtils
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Reproducible science
seed = 12345
np.random.seed(seed)


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


# Function to plot cylinder border
def getCylinderPoints(xc, zc, r):
    xLocOrig1 = np.arange(-r, r + r / 10., r / 10.)
    xLocOrig2 = np.arange(r, -r - r / 10., -r / 10.)
    # Top half of cylinder
    zLoc1 = np.sqrt(-xLocOrig1**2. + r**2.) + zc
    # Bottom half of cylinder
    zLoc2 = -np.sqrt(-xLocOrig2**2. + r**2.) + zc
    # Shift from x = 0 to xc
    xLoc1 = xLocOrig1 + xc * np.ones_like(xLocOrig1)
    xLoc2 = xLocOrig2 + xc * np.ones_like(xLocOrig2)

    topHalf = np.vstack([xLoc1, zLoc1]).T
    topHalf = topHalf[0:-1, :]
    bottomHalf = np.vstack([xLoc2, zLoc2]).T
    bottomHalf = bottomHalf[0:-1, :]

    cylinderPoints = np.vstack([topHalf, bottomHalf])
    cylinderPoints = np.vstack([cylinderPoints, topHalf[0, :]])
    return cylinderPoints


def plot_pseudoSection(
    dc_survey, ax=None, survey_type='dipole-dipole',
    data_type="appConductivity", space_type='half-space',
    clim=None, scale="linear", sameratio=True,
    pcolorOpts={}, data_location=False, dobs=None, dim=2
):
    """
        Read list of 2D tx-rx location and plot a speudo-section of apparent
        resistivity.

        Assumes flat topo for now...

        Input:
        :param SimPEG.EM.Static.DC.SurveyDC.Survey dc_survey: DC survey object
        :param matplotlib.pyplot.axes ax: figure axes on which to plot
        :param str survey_type: Either 'dipole-dipole' | 'pole-dipole' |
            'dipole-pole' | 'pole-pole'
        :param str data_type: Either 'appResistivity' | 'appConductivity' |
            'volt' (potential)
        :param str space_type: Either 'half-space' (default) or 'whole-space'
        :param str scale: Either 'linear' (default) or 'log'

        Output:
        :return  matplotlib.pyplot.figure plot overlayed on image
    """
    import pylab as plt
    from scipy.interpolate import griddata
    # Set depth to 0 for now
    z0 = 0.
    rho = []

    # Use dobs in survey if dobs is None
    if dobs is None:
        if dc_survey.dobs is None:
            raise Exception()
        else:
            dobs = dc_survey.dobs

    rhoApp = DCUtils.apparent_resistivity(
                dc_survey, dobs=dobs,
                survey_type=survey_type,
                space_type=space_type
    )
    midx, midz = DCUtils.source_receiver_midpoints(
                    dc_survey,
                    survey_type=survey_type,
                    dim=dim
    )

    if data_type == 'volt':
        if scale == "linear":
            rho = dobs
        elif scale == "log":
            rho = np.log10(abs(dobs))

    elif data_type == 'appConductivity':
        if scale == "linear":
            rho = 1./rhoApp
        elif scale == "log":
            rho = np.log10(1./rhoApp)

    elif data_type == 'appResistivity':
        if scale == "linear":
            rho = rhoApp
        elif scale == "log":
            rho = np.log10(rhoApp)

    else:
        print()
        raise Exception(
                """data_type must be 'appResistivity' |
                'appConductivity' | 'volt' """
                " not {}".format(data_type)
        )

    # Grid points
    grid_x, grid_z = np.meshgrid(np.linspace(np.min(midx),np.max(midx),int(2.*(np.abs(np.max(midx)-np.min(midx)))+1)),
                              np.linspace(np.min(midz),np.max(midz),int(2.*(np.abs(np.min(midz)-np.max(midz)))+1)),
                               indexing='xy')
    grid_x, grid_z = grid_x.T, grid_z.T

    grid_rho = griddata(np.c_[midx, midz], rho.T, (grid_x, grid_z),
                        method='linear',rescale=True)

    if clim is None:
        vmin, vmax = rho.min(), rho.max()
    else:
        vmin, vmax = clim[0], clim[1]

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(15, 3))

    grid_rho = np.ma.masked_where(np.isnan(grid_rho), grid_rho)
    ph = ax.contourf(
         grid_x[:, 0], grid_z[0, :], grid_rho.T,
        vmin=vmin, vmax=vmax, **pcolorOpts
    )
    ax.set_xlabel('X(m)',fontsize=16)
    ax.set_ylabel('n',fontsize=16)
    ax.tick_params(labelsize=14)
    if scale == "log":
        cbar = plt.colorbar(
               ph, format="$10^{%.1f}$",
               fraction=0.04, orientation="horizontal"
        )
    elif scale == "linear":
        cbar = plt.colorbar(
               ph, format="%.1f",
               fraction=0.04, orientation="horizontal"
        )

    if data_type == 'appConductivity':
        cbar.set_label("App.Cond", size=16)

    elif data_type == 'appResistivity':
        cbar.set_label("Apparent Resistivity (Ohm-m)", size=16)

    elif data_type == 'volt':
        cbar.set_label("Potential (V)", size=16)

    #cmin, cmax = cbar.get_clim()
    #ticks = np.linspace(cmin, cmax, 5)
    #cbar.set_ticks(ticks)
    #cbar.ax.tick_params(labelsize=14)

    cbar.set_ticks(ph.levels)
    cbar.ax.tick_params(labelsize=16)
    # Plot apparent resistivity
    if data_location:
        ax.plot(midx, midz, 'k.', ms=2, alpha=0.8,label='data locations')
        #ax.plot(grid_x, grid_z, 'k.', ms=1, alpha=0.4)
        ax.legend(fontsize=14)
    if sameratio:
        ax.set_aspect('equal', adjustable='box')

    return ax, midx, midz,grid_x,grid_z,rho
