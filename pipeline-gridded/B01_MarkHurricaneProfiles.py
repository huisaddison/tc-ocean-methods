import numpy as np
import pandas as pd
import pickle as pkl

track_dir = './Tracks/'
AL = pd.read_csv(f'{track_dir}/HURDAT_ATLANTIC.csv')
del AL['Unnamed: 0']
EP = pd.read_csv(f'{track_dir}/HURDAT_PACIFIC.csv')
del EP['Unnamed: 0']
WP = pd.read_csv(f'{track_dir}/JTWC_WESTPACIFIC.csv')
del WP['Unnamed: 0']
SH = pd.read_csv(f'{track_dir}/JTWC_SOUTHERNHEMISPHERE.csv')
del SH['Unnamed: 0']
IO = pd.read_csv(f'{track_dir}/JTWC_INDIANOCEAN.csv')
del IO['Unnamed: 0']

Hurricanes = pd.concat([AL, EP, WP, SH, IO])
Hurricanes = Hurricanes[Hurricanes['SEASON'] >= 2007]
Hurricanes = Hurricanes.reset_index(drop=True)
Hurricanes['Timestamp_PD'] = Hurricanes['TIMESTAMP'].apply(
        lambda x: pd.Timestamp(x))
Profiles = pkl.load(open('./Data/ArgoProfileDF.pkl', 'rb'))

near_hur = pd.Series(np.repeat(0, Profiles.shape[0]), dtype=int)

for idx, row in Hurricanes.iterrows():
    if (idx % 1000 == 0):
        print(f'Now processing {idx}')
    dist = (  (Profiles['Latitude']  - row['LAT'])  ** 2
            + (Profiles['Longitude'] - row['LONG']) ** 2).values
    tdiff = (Profiles['Timestamp'] - row['Timestamp_PD'])
    near_hur += pd.Series((
            (tdiff >= pd.Timedelta(days=-30)).values *
            (tdiff <= pd.Timedelta(days=2)).values   *
            (dist <= 64)
        ), dtype=int)

Profiles['NearHurricane'] = near_hur
df_sorted = Profiles.sort_values('Timestamp', ascending=False)
df = df_sorted.drop_duplicates('ProfileID', keep='first')
df = df.reset_index(drop=True)
pkl.dump(df, open('./Data/ArgoProfileDF_NearHur.pkl', 'wb'))

df_no_hur = df[df['NearHurricane'] == 0]
df = df_no_hur.reset_index(drop=True)
pkl.dump(df, open('./Data/ArgoProfileDF_NoHur.pkl', 'wb'))

