import numpy as np
import scipy.sparse as sps
import warnings

from collections import namedtuple
from numba import jit, njit, prange
from scipy.spatial import distance_matrix

from implementation_tools import grid

NULL = object()

# @njit(parallel=True)
def _par_predict_ker(train_X, train_y, X, h):
    n_test = X.shape[0]
    p = X.shape[1]
    test_preds = np.zeros(n_test)
    for idx in prange(n_test):
        k = np.exp(-np.sum((train_X - X[idx, ])**2, axis=1)
                   / (2 * (h**2))) 
        test_preds[idx] = np.dot(k, train_y) / np.sum(k)
    return test_preds

@njit(parallel=True)
def _par_predict_tps(train_X, theta, X, n_test, p):
    preds = np.zeros(n_test)
    for idx in prange(n_test):
        A = np.zeros(p)
        D = np.sum((train_X - X[idx, ]) ** 2, axis=1)
        L = np.log(D)
        L[L == -np.inf] = 0
        E = (1/16) * D * L
        A[:-3] = E
        A[-3] = 1
        A[-2:] = X[idx, ]
        preds[idx] = np.dot(A, theta)
    return preds

class KernelSmoother:
    def __init__(self, h):
        self.train_X = NULL
        self.train_y = NULL
        self.h = h
        return

    def fit(self, train_X, train_y):
        self.train_X = train_X
        self.train_y = train_y

    def predict(self, X):
        return _par_predict_ker(self.train_X, self.train_y, X, self.h)

def _thin_plate_basis(X1, X2=NULL):
    if X2 is NULL:
        X2 = X1
        D = distance_matrix(X2, X1) ** 2
        np.fill_diagonal(D, 1)
    else:
        D = distance_matrix(X2, X1) ** 2
    # It is known that there will be log(0)
    with np.errstate(divide='ignore'):
        L = np.log(D)
        L[L == -np.inf] = 0
    E = (1/16) * (D) * L
    T = np.vstack((np.ones(X2.shape[0]), X2.T))
    return E, T

class ThinPlateSpline:
    def __init__(self, lamb, knots=grid((-8, 8), (-2, 20), 17, 23)):
        self.lamb = lamb
        self.knots = knots
        self.train_X = NULL
        self.train_y = NULL

    def _multiply(self, W, A):
        '''
        Computes the matrix multiplication WA, where W is possibly:

        - NULL (in which we default to the identity)
        - One-dimensional (in which we interpret it as a diagonal matrix)
        - Two-dimensional np.ndarray 
        - Two-dimensional scipy sparse matrix.

        Parameters
        ----------
        W: object or np.ndarray or sps matrix
            Matrix to multiply
        A: np.ndarray
            Matrix to multiply

        Returns
        -------
        np.ndarray
            Result of multiplication.
        '''
        if W is NULL:
            # Interpret as identity
            WA = A
        elif sps.issparse(W):
            WA = W.dot(A)
        # Assume np.ndarray
        elif len(W.shape) == 1:
            # Use array broadcasting
            WA = (W * A.T).T
        else:
            warnings.warn('Dense matrix multiplication.  Use this only if '
                    'matrices are truly dense; otherwise use sparse matrix.')
            WA = np.dot(W, A)
        assert type(WA) is np.ndarray
        return WA

    def fit(self, train_X, train_y, W=NULL):
        '''
        Fit the thin plate spline.

        Parameters
        ----------
        train_X: np.ndarray
            Training design matrix
        train_y: np.ndarray
            Training targets
        W: np.ndarray or sps.sparse_matrix
            Weights for the WLS fit.  Customarily equal to the inverse
            variance.

        Returns
        -------
        None
        '''
        self.train_X = train_X
        self.train_y = train_y
        self.W = W
        # Form basis
        E, T = _thin_plate_basis(self.knots, self.train_X)
        knotE, _ = _thin_plate_basis(self.knots, self.knots)
        self.WE = self._multiply(W, E)
        self.TW = self._multiply(W, T.T).T
        self.A = np.block([
            [np.dot(E.T, self.WE)+self.lamb*knotE,  np.dot(E.T, self.TW.T)],
            [np.dot(T, self.WE),                    np.dot(self.TW, T.T)],
        ])
        self.B = np.vstack((
            self.WE.T,
            self.TW,
        ))
        b = np.dot(self.B, self.train_y)
        self.theta = np.linalg.solve(self.A, b)

    def predict(self, test_X, sd=NULL, S=NULL, k=NULL, diag_H=False):
        '''
        Predict at new locations test_X.

        Parameters
        ----------
        test_X: np.ndarray
            Locations at which to produce predictions.
        sd: bool
            Whether or not to return standard deviations.
        S: np.ndarray or sps.sparse_matrix
            Variance of observations.  Used in the standard deviation
            calculation.

        Returns
        -------
        np.ndarray
            Predictions
        '''
        E, T = _thin_plate_basis(self.knots, test_X)
        C = np.hstack((E, T.T))
        self._test_basis = C
        try:
            preds = np.dot(C, self.theta)
            if sd is NULL:
                ret = preds
            else:
                ret = [None for _ in range(4)]
                ret[0] = preds
                AinvB = np.linalg.solve(self.A, self.B)
                SBtAinv = self._multiply(S, AinvB.T)
                D = np.dot(AinvB, SBtAinv)
                self._cov_theta = D
                if sd == 'diag':
                    var = np.zeros_like(preds)
                    for idx in range(len(preds)):
                        c = C[idx, :]
                        var[idx] = np.dot(c, D.dot(c))
                    sd = np.sqrt(var)
                    ret[1] = sd
                elif sd == 'full':
                    cov = np.linalg.multi_dot((C, D, C.T))
                    var = np.diag(cov)
                    sd = np.sqrt(var)
                    ret[1] = cov
                if k is not NULL:
                    ge = preds - k*sd > 0
                    le = preds + k*sd < 0
                    mask = ge | le
                    ret[2] = mask
                if diag_H:
                    if (test_X.shape[0] != self.train_X.shape[0]):
                        raise ValueError('Evaluation set has different number '
                                'of observations than train set.  Hat matrix '
                                'trick only makes sense when evaluating on '
                                'training set.')
                    Ainv = np.linalg.inv(self.A)
                    Hii = np.zeros(self.train_X.shape[0])
                    for idx in range(self.train_X.shape[0]):
                        Hii[idx] = np.linalg.multi_dot(
                                (C[idx, :], Ainv, self.B[:, idx]))
                    ret[3] = Hii
        except AttributeError as e:
            print(e)
            print("Has tps.train() been called yet?")
        return ret

