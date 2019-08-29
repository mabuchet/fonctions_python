# -*- coding: utf-8 -*-

""" 
Auteur : Marc-Antoine BUCHET

Fonctions utiles à l'ajustement de données par la méthode des moindre carrés.
"""

import numpy as np
from scipy import optimize

################################################################################
# Regression linéaire par moindre carrés :
################################################################################
def reg_lin(x,y,sigma_y):
    """ Calcule les meilleurs paramètres de la droite y=ax+b et leurs incertitudes
    en tenant compte des incertitudes sur y.
    x,y,sigma : listes ou tableaux 1D """
    # version du 28_08_2019
    x=np.array(x)
    y=np.array(y)
    sigma_y=np.array(sigma_y)
    
    w=(1./sigma_y)**2
    Sw=sum(w)
    Sx=sum(w*x)
    Sy=sum(w*y)
    Sxy=sum(w*x*y)
    Sxx=sum(w*x*x)
    
    delta = Sw*Sxx-Sx**2
    a =(Sw*Sxy-Sx*Sy)/delta 
    b = (Sxx*Sy-Sx*Sxy)/delta
    da = np.sqrt(Sw/delta)
    db = np.sqrt(Sxx/delta)
    chisq = sum(w*((y-(a*x+b))**2))
    chisq_red = chisq/(np.size(x)-2.)
    
    return a,da,b,db,chisq,chisq_red

################################################################################
# Moindres carrés :
################################################################################
def leastsq_wrapper(f, xdata, ydata, sigma, p0, ftol=1.49012e-8,
                    xtol=1.49012e-8, maxfev=0,full_output=0):
    """
    Wrapper for leastsq drawn from curve_fit but adapted a bit.
    Use non-linear least squares to fit a function, f, to data.
    Assumes ``ydata = f(xdata, *params) + eps``
    
    Parameters
    ----------
    f : callable
        The model function, f(x, ...). It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    xdata : An M-length sequence or an (k,M)-shaped array
        for functions with k predictors.
        The independent variable where the data is measured.
    ydata : M-length sequence
        The dependent data --- nominally f(xdata, ...)
    sigma : None or M-length sequence, these values are used as weights in the
        least-squares problem. `sigma` should describe one standard deviation
        errors of the input data points. The estimated covariance in `pcov` is
        based on these values.
    p0 : scalar, or N-length sequence
        Initial guess for the parameters.
        
    full_output : bool
        non-zero to return all optional outputs.
    ftol : float
        Relative error desired in the sum of squares.
    xtol : float
        Relative error desired in the approximate solution.
    maxfev : int
        The maximum number of calls to the function. If zero, then 100*(N+1) is
        the maximum where N is the number of elements in x0.
        
    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the squared error
        of ``f(xdata, *popt) - ydata`` is minimized
    pcov : 2d array
        The estimated covariance of popt. The diagonals provide the variance
        of the parameter estimate. To compute one standard deviation errors
        on the parameters use ``perr = np.sqrt(np.diag(pcov))``.
    chi_sq : scalar
        Value of the chi-square fonction evaluated with popt.
    infodict : dict
        a dictionary of optional outputs with the key s:
        
        ``nfev``
            The number of function calls
        ``fvec``
            The function evaluated at the output
        ``fjac``
            A permutation of the R matrix of a QR
            factorization of the final approximate
            Jacobian matrix, stored column wise.
            Together with ipvt, the covariance of the
            estimate can be approximated.
        ``ipvt``
            An integer array of length N which defines
            a permutation matrix, p, such that
            fjac*p = q*r, where r is upper triangular
            with diagonal elements of nonincreasing
            magnitude. Column j of p is column ipvt(j)
            of the identity matrix.
        ``qtf``
            The vector (transpose(q) * fvec).
            
    mesg : str
        A string message giving information about the cause of failure.
    ier : int
        An integer flag. If it is equal to 1, 2, 3 or 4, the solution was
        found. Otherwise, the solution was not found. In either case, the
        optional output variable 'mesg' gives more information.
    
    See Also
    --------
    leastsq
    
    Notes
    -----
    The algorithm uses the Levenberg-Marquardt algorithm through `leastsq`.
    Additional keyword arguments are passed directly to that algorithm.
    """
    # version du 28_08_2019
    def func(params, xdata, ydata, function, weights):
        return weights * (function(xdata, *params) - ydata)
    
    # Check input arguments
    if np.isscalar(p0):
        p0 = np.array([p0])
    
    xdata = np.array(xdata).flatten()
    ydata = np.array(ydata).flatten()
    
    
    args = (xdata, ydata, f,1.0 / np.asarray(sigma))
    res = optimize.leastsq(func, p0, args=args, full_output=1, ftol=ftol, xtol=xtol,
                           maxfev=maxfev)
    (popt, pcov, infodict, errmsg, ier) = res
    
    if ier not in [1, 2, 3, 4]:
        msg = "Optimal parameters not found: " + errmsg
        raise RuntimeError(msg)
    chisq = (np.asarray(func(popt, *args))**2).sum()
    if pcov is None:
        # indeterminate covariance
        pcov = np.zeros((len(popt), len(popt)), dtype=float)
        pcov.fill(inf)
    perr = np.sqrt(np.diag(pcov))
    chisq_red=chisq/(len(xdata)-len(p0))
    if full_output :
        return popt, perr, chisq, chisq_red, (pcov, infodict, errmsg, ier)
    else : return popt, perr, chisq, chisq_red