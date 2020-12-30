import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
# from mpl_toolkits.basemap import Basemap
from functools import partial
from tools import (
        haversine,
        pd_matlab_days,
        matlab_to_datetime,
        datetime_to_matlab,
        conform_lons,
        geodesic_signed_angle,
        flat_earth_angle,
        flat_earth_signed_angle,
        date_to_yearday,
        create_ArgoProfileID,
        induced_angle,
        angle_sign,
        )
NULL = object()
matlab_to_datetime_v = np.vectorize(matlab_to_datetime)

class Processor:
    def __init__(self, hurricane_df, f):
        self.hurricane_df = hurricane_df
        self.f = f
        self.dates = hurricane_df['DATE'].apply(
            lambda x: pd.to_datetime(x))
        self.cycle_nums = np.array(f['profCycleNumberAggr']).flatten()
        self.float_ids = np.array(f['profFloatIDAggr']).flatten()
        self.argo_dates = np.array(f['profJulDayAggr']).flatten()
        self.argo_datetimes = matlab_to_datetime_v(self.argo_dates)
        self.argo_lons = conform_lons(
                np.array(f['profLongAggr'])).flatten()
        self.argo_lats = np.array(f['profLatAggr']).flatten()
        self.profile_ids = np.array(list(
            create_ArgoProfileID(f, c)
            for f, c in zip(self.float_ids, self.cycle_nums)))
        self.after_floats = set()

    def _bounding_box_floats(self,
            time,
            time_diff_1,
            time_diff_2,
            lon,
            lon_diff_1,
            lon_diff_2,
            lat,
            lat_diff_1,
            lat_diff_2,
            ):
        mn_time = pd_matlab_days(time +
                pd.Timedelta(days=time_diff_1))
        mx_time = pd_matlab_days(time +
                pd.Timedelta(days=time_diff_2))
        mn_lon = lon + lon_diff_1
        mx_lon = lon + lon_diff_2
        mn_lat = lat + lat_diff_1
        mx_lat = lat + lat_diff_2

        mn_date_mask = self.argo_dates > mn_time
        mx_date_mask = self.argo_dates < mx_time
        mn_lon_mask = self.argo_lons > mn_lon
        mx_lon_mask = self.argo_lons < mx_lon
        mn_lat_mask = self.argo_lats > mn_lat
        mx_lat_mask = self.argo_lats < mx_lat
        mask = (mn_date_mask
                * mx_date_mask
                * mn_lon_mask
                * mx_lon_mask
                * mn_lat_mask
                * mx_lat_mask)
        # return list(np.nonzero(mask.flatten())[0])
        return mask

    def _generate_before_float_candidates(self):
        counts = np.zeros_like(self.argo_dates)
        for idx, row in self.hurricane_df.iterrows():
            counts += self._bounding_box_floats(
                time        = pd.Timestamp(row['TIMESTAMP']),
                time_diff_1 = -12,
                time_diff_2 = -2,
                lon         = row['LONG'],
                lon_diff_1  = -10,
                lon_diff_2  = +10,
                lat         = row['LAT'],
                lat_diff_1  = -10,
                lat_diff_2  = +10,
                )
        return self.profile_ids[np.nonzero(counts)[0]]
    
    def _get_lat(self, pid):
        return self.argo_lats[self.profile_ids == pid][0]

    def _get_lon(self, pid):
        return self.argo_lons[self.profile_ids == pid][0]

    def _get_t(self, pid):
        return self.argo_datetimes[self.profile_ids == pid][0]

    def _project_float_onto_track(self, argo_profile_id):
        track_len = self.hurricane_df.shape[0]
        argo_lon = self.argo_lons[self.profile_ids == argo_profile_id][0]
        argo_lat = self.argo_lats[self.profile_ids == argo_profile_id][0]
        argo_t   = matlab_to_datetime(
                self.argo_dates[self.profile_ids == argo_profile_id][0])
        df = self.hurricane_df.assign(
                DistToFloat = lambda r:
                np.sqrt((r.LONG - argo_lon) ** 2 + (r.LAT - argo_lat) ** 2))
        '''
        Establish pair of track coordinates between which the float lies.
        '''
        loc1 = np.argmin(df['DistToFloat'].values)
        if (loc1 == 0) or (loc1 == track_len - 1):
            return None
        loc2 = (loc1 - 1) if (
                df['DistToFloat'].values[loc1-1] < df['DistToFloat'].values[loc1+1]
                ) else (loc1 + 1)
        wind     = df['WIND'].values[loc1]
        hurclass = df['CLASS'].values[loc1]
        lat1, lon1 = df['LAT'].values[loc1], df['LONG'].values[loc1]
        lat2, lon2 = df['LAT'].values[loc2], df['LONG'].values[loc2]
        theta0 = flat_earth_angle(lat1, lon1, lat2, lon2, argo_lat, argo_lon)
        in_cone = True if (theta0 > 90) else False
        along_track_tot = np.sqrt((lat2 - lat1) ** 2 + (lon2 - lon1) ** 2)
        if in_cone is False: # Lies between loc1 and loc2
            # Ensure proper ordering
            loc1, loc2 = (loc1, loc2) if (loc1 < loc2) else (loc2, loc1)
            lat1, lon1 = df['LAT'].values[loc1], df['LONG'].values[loc1]
            lat2, lon2 = df['LAT'].values[loc2], df['LONG'].values[loc2]
            dist1 = df['DistToFloat'].values[loc1]
            theta1 = flat_earth_angle(
                    lat1, lon1,
                    lat2, lon2,
                    argo_lat, argo_lon)
            # projecting location of argo float onto the track
            # along_track = np.sqrt((dist1 ** 2) / (1 + np.tan(theta1) ** 2))
            along_track = np.abs(dist1 * np.cos(np.deg2rad(theta1)))
            alpha = along_track / along_track_tot # \in [0, 1]

            t1 = pd.Timestamp(df['TIMESTAMP'].values[loc1])
            t2 = pd.Timestamp(df['TIMESTAMP'].values[loc2])
            proj_lat = lat1 + alpha * (lat2 - lat1)
            proj_lon = lon1 + alpha * (lon2 - lon1)
            proj_t   = t1   + alpha * (t2   -   t1)
        else:
            # Establish which float to project to
            alpha, proj_loc, loc1, loc2 = (
                    0, loc1, loc1, loc2) if (loc1 < loc2) else (
                    1, loc1, loc2, loc1)
            along_track = along_track_tot * alpha
            proj_lat = df['LAT'].values[proj_loc]
            proj_lon = df['LONG'].values[proj_loc]
            proj_t   = pd.Timestamp(df['TIMESTAMP'].values[proj_loc])
        angle    = induced_angle(argo_lat, argo_lon, proj_lat, proj_lon)
        sign     = angle_sign(proj_lat, proj_lon, lat2, lon2, argo_lat, argo_lon)
        if ((angle > 8)
                or (argo_t - proj_t < pd.Timedelta(days=-12))
                or (argo_t - proj_t > pd.Timedelta(days=-2))
            ):
            return None
        else:
            return {
                'before_pid':               argo_profile_id,
                'argo_t':                   argo_t,
                'argo_lat':                 argo_lat,
                'argo_lon':                 argo_lon,
                'along_track':              along_track,
                'along_track_tot':          along_track_tot,
                'alpha':                    alpha,
                'proj_lat':                 proj_lat,
                'proj_lon':                 proj_lon,
                'proj_t':                   proj_t,
                'loc1':                     loc1,
                'loc2':                     loc2,
                'angle':                    angle,
                'sign':                     sign,
                'wind':                     wind,
                'class':                    hurclass,
            }

    def generate_before_floats(self):
        candidates = self._generate_before_float_candidates()
        row_list = []
        for pid in candidates:
            proj = self._project_float_onto_track(pid)
            if proj is not None:
                row_list.append(proj)
        if len(row_list) > 0:
            self.float_df = pd.DataFrame(row_list)
        else:
            self.float_df = pd.DataFrame(columns=[
                'before_pid',
                'argo_t',
                'argo_lat',
                'argo_lon',
                'along_track',
                'along_track_tot',
                'alpha',
                'proj_lat',
                'proj_lon',
                'proj_t',
                'loc1',
                'loc2',
                'angle',
                'sign',
                'wind',
                'class',
            ])

    def _generate_after_floats(self, row, td1=-2, td2=3):
        mask = self._bounding_box_floats(
            time        = row['proj_t'],
            time_diff_1 = td1,
            time_diff_2 = td2,
            lon         = row['argo_lon'],
            lon_diff_1  = -0.5,
            lon_diff_2  = +0.5,
            lat         = row['argo_lat'],
            lat_diff_1  = -0.5,
            lat_diff_2  = +0.5,
            )
        candidate_pids = self.profile_ids[np.nonzero(mask)[0]]
        after_floats = set(cpid for cpid in candidate_pids if
                (induced_angle(
                    row['argo_lat'],     row['argo_lon'],
                    self._get_lat(cpid), self._get_lon(cpid)) < 0.2)
        )
        # after_floats.discard(row['before_pid'])
        self.after_floats.update(after_floats)
        return after_floats

    def add_after_floats(self):
        self.float_df = self.float_df.assign(
                AfterForced=self.float_df.apply(lambda r:
                    self._generate_after_floats(r, -2,  3), axis=1),
                AfterRestore=self.float_df.apply(lambda r:
                    self._generate_after_floats(r, 3, 20), axis=1),
                )
        self.float_df = self.float_df.assign(
                AfterCombined=self.float_df.apply(lambda r:
                    r['AfterForced'].union(r['AfterRestore']), axis=1)
                )
        self.float_df = self.float_df.assign(
                AfterForcedNum=self.float_df['AfterForced'].apply(len),
                AfterRestoreNum=self.float_df['AfterRestore'].apply(len),
                AfterCombinedNum=self.float_df['AfterCombined'].apply(len),
                )

    def _get_float_pairs(self, stage='Forced'):
        pairs = []
        df = self.float_df[self.float_df[f'After{stage}Num'] > 0][[
                'before_pid',
                'argo_t', # 'Before' datetime
                f'After{stage}',
                'angle',
                'sign',
                'wind',
                'proj_t',
                ]]
        for idx, row in df.iterrows():
            pairs.extend([(
                row['before_pid'],
                row['argo_t'], # 'Before' datetime
                x,
                self._get_t(x), # 'After' datetime
                row['angle'],
                row['sign'],
                row['wind'],
                row['proj_t'],
                )
                for x in row[f'After{stage}']])
        return pairs

    def create_pair_df(self):
        stage = 'Combined'
        df = self.float_df[self.float_df[f'After{stage}Num'] > 0][[
                'before_pid',
                'argo_t', # 'Before' datetime
                f'After{stage}',
                f'After{stage}Num',
                'angle',
                'sign',
                'wind',
                'proj_t',
                ]]
        df_lst = []
        for idx, row in df.iterrows():
            n = row[f'After{stage}Num']
            after_floats = row[f'After{stage}']
            df_lst.append(pd.DataFrame({
                'before_pid':   [row['before_pid'] for _ in range(n)],
                'before_t':     [row['argo_t'] for _ in range(n)],
                'after_pid':    [pid for pid in after_floats],
                'after_t':      [self._get_t(pid) for pid in after_floats],
                'angle':        [row['angle'] for _ in range(n)],
                'wind':         [row['wind'] for _ in range(n)],
                'proj_t':       [row['proj_t'] for _ in range(n)],
                'sign':         [row['sign'] for _ in range(n)],
                }))

        if len(df_lst) > 0:
            self.float_pairs_df = pd.concat(df_lst
                    ).sort_values('before_t', ascending=False
                    ).drop_duplicates('after_pid' # Keep only pair w/ newest before
                    ).reset_index(drop=True)
        else:
            self.float_pairs_df = None
        return self.float_pairs_df

    def _get_temperature_diffs(self, temp_dict, stage='Forced'):
        def minus(a, b):
            if a is None or b is None:
                return None
            else:
                return a - b
        pairs = self._get_float_pairs(stage=stage)
        td = pkl.load(open(temp_dict, 'rb'))
        diffs = np.fromiter((minus(td.get(after, None), td.get(before, None)) for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs), dtype=float)
        angles = np.fromiter((- sign * angle for # Negative sign to flip left right
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs), dtype=float)
        winds  = np.fromiter((wind for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs), dtype=float)
        before= np.fromiter((before for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs), dtype=np.dtype(('U', 16)))
        after = np.fromiter((after for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs), dtype=np.dtype(('U', 16)))
        before_times_lst = [b_t for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs]
        after_times_lst = [a_t for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs]
        hurricane_times_lst = [p_t for
            (before, b_t, after, a_t, angle, sign, wind, p_t) in pairs]
        before_times = np.array(before_times_lst)
        after_times = np.array(after_times_lst)
        hurricane_times = np.array(hurricane_times_lst)
        return (
            before[~np.isnan(diffs)],
            after[~np.isnan(diffs)],
            diffs[~np.isnan(diffs)],
            angles[~np.isnan(diffs)],
            winds[~np.isnan(diffs)],
            before_times[~np.isnan(diffs)],
            after_times[~np.isnan(diffs)],
            hurricane_times[~np.isnan(diffs)],
            )

    def plot_pair_differences(self, temp_dict, stage='Forced', bins=20, min_wind=NULL):
        diffs, angles, winds = self._get_temperature_diffs(stage=stage,
                temp_dict=temp_dict)
        if min_wind is not NULL:
            plt.hist(diffs[winds >= min_wind], bins=bins)
        else:
            plt.hist(diffs, bins=bins)

    def plot(self, basin='north_atlantic', after_floats=False, float_labels=False):
        if basin=='north_atlantic':
            bm_params = {
                'llcrnrlon':-100.,
                'llcrnrlat':-5.,
                'urcrnrlon':20.,
                'urcrnrlat':50.,
                'projection':'lcc',
                'lat_1':20.,
                'lat_2':40.,
                'lon_0':-60.,
                'resolution' :'l',
                'area_thresh':1000.,
            }
            dp_params = {
                'circles' : np.arange(10, 70, 20),
                'labels' : [1, 1, 0, 0],
            }
            dm_params = {
                'meridians' : np.arange(-100,0,20),
                'labels' : [0,0,0,1],
            }
            dms_params = {
                'lon' : -90,
                'lat' : 3,
                'lon0' : 60,
                'lat0' : 30,
                'length' : 2000,
                'barstyle' : 'fancy',
            }
        elif basin=='east_pacific':
            bm_params = {
                'llcrnrlon':-160.,
                'llcrnrlat':-0.,
                'urcrnrlon':-40.,
                'urcrnrlat':35.,
                'projection':'lcc',
                'lat_1':20.,
                'lat_2':30.,
                'lon_0':-150.,
                'resolution' :'l',
                'area_thresh':1000.,
            }
            dp_params = {
                'circles' : np.arange(00, 40, 20),
                'labels' : [1, 1, 0, 0],
            }
            dm_params = {
                'meridians' : np.arange(-160,-40,20),
                'labels' : [0,0,0,1],
            }
            dms_params = {
                'lon' : -150,
                'lat' : 5,
                'lon0' : -120,
                'lat0' : 17,
                'length' : 2000,
                'barstyle' : 'fancy',
            }
        elif basin=='west_pacific':
            bm_params = {
                'llcrnrlon':110.,
                'llcrnrlat':0.,
                'urcrnrlon':200.,
                'urcrnrlat':35.,
                'projection':'lcc',
                'lat_1':20.,
                'lat_2':40.,
                'lon_0':140.,
                'resolution' :'l',
                'area_thresh':1000.,
            }
            dp_params = {
                'circles' : np.arange(00, 40, 10),
                'labels' : [1, 1, 0, 0],
            }
            dm_params = {
                'meridians' : np.arange(110,200,20),
                'labels' : [0,0,0,1],
            }
            dms_params = {
                'lon' : 172,
                'lat' : 3,
                'lon0' : 140,
                'lat0' : 30,
                'length' : 2000,
                'barstyle' : 'fancy',
            }
        else:
            raise Exception('Invalid basin specified.')
        m = Basemap(**bm_params)
        m.drawcoastlines()
        # m.drawcountries() 
        m.drawmapboundary(fill_color='#99ffff')
        m.fillcontinents(color='#cc9966',lake_color='#99ffff', zorder=1)
        m.drawparallels(**dp_params)
        m.drawmeridians(**dm_params)
        m.drawmapscale(**dms_params)
        x, y = m(np.array(self.hurricane_df['LONG']),
                 np.array(self.hurricane_df['LAT']))
        plt.plot(x, y,
                color='r', linewidth=0.5)
        plt.scatter(x, y, color='r', marker='*', s=0.5)
        for idx, r in self.float_df.iterrows():
            xa, ya = m(r['argo_lon'], r['argo_lat'])
            xp, yp = m(r['proj_lon'], r['proj_lat'])
            _ = plt.scatter(xa, ya, color='k', s=0.5)
            _ = plt.plot((xp, xa), (yp, ya), color='k', linewidth=0.5, alpha=0.5)
            if float_labels:
                plt.annotate(
                    f"{r['before_pid']}",
                    (xa, ya), color='k', alpha=0.5)
        if after_floats:
            for apid in self.after_floats:
                x, y = m(self._get_lon(apid), self._get_lat(apid))
                plt.scatter(x, y, color='tab:orange', s=2.5)
            if float_labels:
                plt.annotate(
                    f"{apid}",
                    (x, y), color='k', alpha=0.5)


class BackgroundProcessor(Processor):
    def __init__(self, hurricane_df, df, year=2016):
        hurricane_df = hurricane_df.reset_index(drop=True)
        hurricane_df['TIMESTAMP'] = hurricane_df['TIMESTAMP'].apply(
                lambda x: pd.to_datetime(x).replace(year=year))
        hurricane_df['DATE'] = hurricane_df['DATE'].apply(
                lambda x: pd.to_datetime(x).replace(year=year))
        self.hurricane_df = hurricane_df
        self.argo_df = df
        self.dates = hurricane_df['DATE']
        self.cycle_nums = self.argo_df['CycleNum'].values
        self.float_ids =  self.argo_df['FloatID'].values
        self.argo_dates = self.argo_df['ArgoDate'].values
        self.argo_datetimes = self.argo_df['Timestamp'].values
        self.argo_lons =  self.argo_df['Longitude'].values
        self.argo_lats =  self.argo_df['Latitude'].values
        self.profile_ids = self.argo_df['ProfileID'].values 
        self.after_floats = set()
