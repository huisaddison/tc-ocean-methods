import typing
import numpy as np
import pandas as pd
import warnings
from math import radians, cos, sin, asin, sqrt
from datetime import date as ddate, datetime, timedelta

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points (lon1, lat1), (lon2,
    lat2), where longitude and latitude are specified in degrees.

    Reference: https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def pd_matlab_days(pd_dt):
    """
    Converts Pandas datetime object to Matlab days since
    year zero.
    """
    return datetime_to_matlab(pd_dt.to_pydatetime())

def datetime_to_matlab(dt):
    """
    Converts Datetime object to Matlab days since year zero.
    """
    return ((dt - datetime(1, 1, 1)).total_seconds() /
        timedelta(days=1).total_seconds()) + 367

def matlab_to_datetime(md):
    """
    Converts Matlab days since year zero to Python
    datetime object.
    """
    return datetime(1, 1, 1) + timedelta(md) - timedelta(367)

def matlab_to_datetime_v(md_array):
    return np.vectorize(matlab_to_datetime)

def conform_ghrsst_time(ds):
    return datetime(1981, 1, 1) + timedelta(seconds=ds)

def date_to_yearday(date: tuple) -> int:
    return ddate(*date).timetuple().tm_yday

def conform_lons(lons):
    """
    Maps longitude from [20, 380] to [-180, 180].

    For use with Argo data, which maps longitude onto [20, 380].
    """
    to_fix = lons > 180
    fixed = lons - 360
    return (1 - to_fix) * lons + to_fix * fixed

def lift(lat: float, lon: float) -> np.ndarray:
    """Lifts the lat/lon pair into 3D coordinates."""
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    return np.array([
        np.cos(lat) * np.sin(lon),
        np.cos(lat) * np.cos(lon),
        np.sin(lat)])

def lifted_angle(u: np.ndarray, v: np.ndarray) -> float:
    """Returns the angle between two vectors."""
    inprod = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return np.rad2deg(np.arccos(inprod))

def induced_angle(lat1, lon1, lat2, lon2):
    p1 = lift(lat1, lon1)
    p2 = lift(lat2, lon2)
    return lifted_angle(p1, p2)

def great_circle_angle(
        lat0: float, lon0: float,
        lat1: float, lon1: float,
        lat2: float, lon2: float) -> float:
    """
    Computes the angle between location0 -> location 1
    and location 0 -> location 2 along the geodesic.

    Reference: https://math.stackexchange.com/questions/2809181/how-to-calculate-angle-between-two-directions-on-sphere
    """
    loc0 = lift(lat0, lon0)
    loc1 = lift(lat1, lon1)
    loc2 = lift(lat2, lon2)
    u = np.cross(loc0, loc1)
    v = np.cross(loc0, loc2)
    return lifted_angle(u, v)
    
def geodesic_signed_angle(
        lat0: float, lon0: float,
        lat1: float, lon1: float,
        lat2: float, lon2: float) -> float:
    """
    Computes the angle between the vectors u and v, where:
        u = (lon1 - lon0, lat1 - lat0)
        v = (lon2 - lon0, lat2 - lat0)

    This computation uses the great circle approximation, i.e., the Earth is 
    assumed to be spherical.

    CAUTION: The sign of the angle may be implemented incorrectly and in 
    testing is incorrect half the time.  The magnitudes of the angles are 
    implemented correctly.  We recommend using `flat_earth_signed_angle()`,
    as the discrepancy at the ~1000km scale between the two approximations
    are minimal.
    """
    loc0 = lift(lat0, lon0)
    loc1 = lift(lat1, lon1)
    loc2 = lift(lat2, lon2)
    u = np.cross(loc0, loc1)
    v = np.cross(loc0, loc2)
    cr = np.cross(u, v)
    warnings.warn("Sign of angle may be implemented incorrectly.  Use with caution.")
    sign = -np.sign(np.dot(cr, loc0))
    raw_angle = lifted_angle(u, v)
    signed_angle = raw_angle * sign
    return signed_angle

def flat_earth_angle(
        lat0: float, lon0: float,
        lat1: float, lon1: float,
        lat2: float, lon2: float) -> float:
    u = np.array([lat1 - lat0, lon1 - lon0, 0])
    v = np.array([lat2 - lat0, lon2 - lon0, 0])
    inprod = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    angle =  np.rad2deg(np.arccos(inprod))
    return angle

def flat_earth_signed_angle(
        lat0: float, lon0: float,
        lat1: float, lon1: float,
        lat2: float, lon2: float) -> float:
    """
    Computes the angle between the vectors u and v, where:
        u = (lon1 - lon0, lat1 - lat0)
        v = (lon2 - lon0, lat2 - lat0)

    This computation assumes the Earth to be flat.  The errors due to this 
    approximation are negligible at the ~1000km scale.
    """
    u = np.array([lat1 - lat0, lon1 - lon0, 0])
    v = np.array([lat2 - lat0, lon2 - lon0, 0])
    cr = np.cross(u, v)
    inprod = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    angle =  np.rad2deg(np.arccos(inprod))
    sign = -np.sign(cr[2])
    signed_angle = angle * sign
    return signed_angle

def angle_sign(
        lat0: float, lon0: float,
        lat1: float, lon1: float,
        lat2: float, lon2: float) -> float:
    """
    Computes the angle between the vectors u and v, where:
        u = (lon1 - lon0, lat1 - lat0)
        v = (lon2 - lon0, lat2 - lat0)

    This computation assumes the Earth to be flat.  The errors due to this 
    approximation are negligible at the ~1000km scale.
    """
    u = np.array([lat1 - lat0, lon1 - lon0, 0])
    v = np.array([lat2 - lat0, lon2 - lon0, 0])
    cr = np.cross(u, v)
    sign = -np.sign(cr[2])
    return sign

def create_ArgoProfileID(fid, cn, delta=0):
    return f'{int(fid):010}_{int(cn)+delta:05}'

def replace(infile: str, outfile: str, replacement_dict: dict):
    with open(infile) as f:
        newtext = f.read()
        for k, v in replacement_dict.items():
            newtext = newtext.replace(k, v)
    with open(outfile, 'w') as f:
        f.write(newtext)

def covariance_matrix(df: pd.DataFrame, df_param: pd.DataFrame,
        depth_idx: int):
    '''
    Covariance matrix for the subset of pairs in df.

    All pairs must be associated with the same before float.  If
    df only has one row, then a single value equal to the variance
    of that pair is returned.

    Parameters
    ----------
    df: pd.DataFrame
        Subset of pairs for which to calculate a covariance matrix.
    df_param: pd.DataFrame
        Parameters from the Gaussian Process fit, indexed by a unique
        before_pid.
    depth_idx: int
        Depth index at which to compute covariance.

    Returns
    -------
    float or np.ndarray[float]
    '''
    if len(df['before_pid'].unique()) != 1:
        raise ValueError('Must have exactly one before_pid')
    n_pairs = df.shape[0]
    # Read values
    before_pid = df['before_pid'].values[0]
    r = df_param.loc[before_pid]
    try:
        theta_t = r['theta_t'][depth_idx]
        phi = r['phi'][depth_idx]
        sigma = r['sigma'][depth_idx]
    except TypeError:
        return np.nan
    # Construct covariance matrix for y
    before_t = df['before_t'].values[0]
    t = np.append(before_t, df['after_t'].values)
    t_diff = (t - t.reshape(-1, 1)) / pd.to_timedelta(1, unit='D')
    d = np.abs(t_diff / theta_t)
    K = phi * np.exp(-d)
    np.fill_diagonal(K, phi + sigma**2)
    # Construct covariance matrix for dy
    A = np.hstack((
        -np.ones(n_pairs).reshape(-1, 1),
        np.identity(n_pairs)))
    return np.linalg.multi_dot((A, K, A.T))
