from scipy.constants import mu_0
import numpy as np

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


def omega(frequency):
    """
    angular frequency
    """
    return 2*np.pi*frequency


def appres_phase_from_data(survey):
    """
    Compute apparent resistivity and phase given impedances
    (real and imaginary components)
    and the frequency.
    """
    data = survey.dobs
    frequency = survey.frequency
    # data are arranged (Zxy_real, Zxy_imag) for each frequency
    Zxy_real = data.reshape((survey.nFreq, 2))[:, 0]
    Zxy_imag = data.reshape((survey.nFreq, 2))[:, 1]
    Zxy = Zxy_real+1j*Zxy_imag

    # compute apparent resistivity and phase from complex impedance
    app_res = abs(Zxy)**2 / (mu_0*omega(frequency))
    phase = np.rad2deg(np.arctan(Zxy_imag / Zxy_real))

    return app_res, phase


def weighted_avg_and_var(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2., weights=weights)
    return (average, variance)
