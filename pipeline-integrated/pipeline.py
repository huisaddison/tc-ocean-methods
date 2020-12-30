from datetime import datetime
import os
import subprocess

from tools import replace

def matlab_run(script: str, replacement_dict: dict = {}):
    timestamp = datetime.now()
    f_out = f"temp_matlab_{datetime.now().strftime('%Y%m%d_%H%M%S')}.m"
    replace(script, f_out, replacement_dict)
    proc=subprocess.run([
        'matlab',
        '-nodisplay',
        '-nosplash',
        '-nodesktop',
        f'-r "run(\'{f_out}\');exit;"'
        ])
    if proc.returncode != 0:
        raise RuntimeError(f'Subprocess {f_out} exited with non-zero return '
                'status.')
    # cleanup
    os.remove(f_out)
    return 

## A01
YEARS_LIST = [
        '2007_2010',
        '2011_2014',
        '2015_2016',
        '2017_2018',
]

for YEARS in YEARS_LIST:
    matlab_run('A01_pchipIntegration.m', {
        '<PY:YEARS>': f'{YEARS}',
        '<PY:DATA_LOC>': '$HOME/',
        '<PY:GRID_LOWER>': '10',
        '<PY:GRID_UPPER>': '200',
    })

## A02
matlab_run('A02_concatenateArrays.m')

'''
## A03
matlab_run('A03_createDataMask.m', {
    '<PY:WINDOW_SIZE>': '8',
})
'''

## A04
matlab_run('A04_filterUsingMasks.m', {
    '<PY:WINDOW_SIZE>': '8',
})

## A05: Hurricane profiles
matlab_run('A05_splitHurricaneProfiles.m', {
    '<PY:GRID_TEMP_FN>': './Data/gridTempProfHurricane_',
    '<PY:MASK_VALUE>': '0',
    '<PY:MASK_NAME>': 'NoHurMask.csv',
    '<PY:WINDOW_SIZE>': '8',
})

## A05: Non-hurricane profiles
matlab_run('A05_splitHurricaneProfiles.m', {
    '<PY:GRID_TEMP_FN>': './Data/gridTempProfNonHurricane_',
    '<PY:MASK_VALUE>': '1',
    '<PY:MASK_NAME>': 'NoHurMask.csv',
    '<PY:WINDOW_SIZE>': '8',
}

## A06
matlab_run('A06_estimateMeanField.m', {
    '<PY:START_YEAR>': '2007',
    '<PY:END_YEAR+1>': '2019',
    '<PY:WINDOW_SIZE>': '8',
})

## Plotting
matlab_run('A13_plotMeanField.m', {
    '<PY:WINDOW_SIZE>': '8',
})

## A07: subtractMean for hurricanes
matlab_run('A07_subtractMean.m', {
    '<PY:GRID_TEMP_FN>':        './Data/gridTempProfFiltered_',
    '<PY:RES_TEMP_FN>':         './Data/gridTempRes_',
    '<PY:WINDOW_SIZE>': '8',
})

## A07: subtractMean for non-hurricanes
matlab_run('A07_subtractMean.m', {
    '<PY:GRID_TEMP_FN>':        './Data/gridTempProfNonHurricane_',
    '<PY:RES_TEMP_FN>':         './Data/gridTempResNonHurricane_',
    '<PY:WINDOW_SIZE>': '8',
})

## A08
matlab_run('A08_divideDataToMonths.m', {
    '<PY:START_YEAR>': '2007',
    '<PY:END_YEAR>': '2018',
    '<PY:WINDOW_SIZE>': '8',
})

## A09
matlab_run('A09_extendedData.m', {
    '<PY:START_YEAR>': '2007',
    '<PY:END_YEAR>': '2018',
    '<PY:WINDOW_SIZE>': '8',
})

## A10: WP
matlab_run('A10_filterLocalMLESpaceTime.m', {
    '<PY:START_YEAR>': '2007',
    '<PY:END_YEAR>': '2018',
    '<PY:WINDOW_SIZE>': '8',
    '<PY:CENTER_MONTH>': '9',
    '<PY:OCEAN_BASIN>': '_WestPacific',
})

## A10: NA
matlab_run('A10_filterLocalMLESpaceTime.m', {
    '<PY:START_YEAR>': '2007',
    '<PY:END_YEAR>': '2018',
    '<PY:WINDOW_SIZE>': '8',
    '<PY:CENTER_MONTH>': '9',
    '<PY:OCEAN_BASIN>': '_NorthAtlantic',
})

# A11
OB = [
    # ('_NorthAtlantic', 'meshgrid(0.5:70.5,261.5:360.5)'),
    # ('_WestPacific',   'meshgrid(0.5:65.5,105.5:187.5)'),
    ('_AllBasins',     'meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360))'),
]

start_year = '2007'
end_year = '2018'
window_size = '8'
window_size_gp = '8'
center_month = '9'
fn = 'Results/localMLESpaceTime_Depth_{d:03d}_{ws}_20_{window_size_gp}_{cm:02d}_{sy}_{ey}{ob}.mat'
for ob, ob_mesh in OB:
    if not os.path.exists(fn.format(
        ws=window_size,
        cm=int(center_month),
        sy=start_year,
        ey=end_year,
        ob=ob)):
        print(ob)
        matlab_run('A11_localMLESpaceTime.m', {
            '<PY:START_YEAR>': start_year,
            '<PY:END_YEAR>': end_year,
            '<PY:WINDOW_SIZE>': window_size,
            '<PY:WINDOW_SIZE_GP>': window_size_gp,
            '<PY:CENTER_MONTH>': center_month,
            '<PY:OCEAN_BASIN>': ob,
            '<PY:OB_MESHGRID>': ob_mesh,
        })

# A12 - NA
OB = [
    # ('_NorthAtlantic', 'meshgrid(0.5:70.5,261.5:360.5)'),
    # ('_WestPacific',   'meshgrid(0.5:65.5,105.5:187.5)'),
    ('_AllBasins',     'meshgrid(linspace(-89.5,89.5,180),linspace(20.5,379.5,360))'),
]

start_year = '2007'
end_year = '2018'
window_size = '8'
window_size_gp = '8'
center_month = '9'
fn = 'Results/localMLESpaceTime_Depth_{d:03d}_{ws}_20_{window_size_gp}_{cm:02d}_{sy}_{ey}{ob}.mat'
for ob, ob_mesh in OB:
    print(ob)
    matlab_run('A12_fitLocalMLESpaceTime.m', {
        '<PY:GRID_TEMP_FN>': './Data/gridTempProfHurricane_',
        '<PY:START_YEAR>': start_year,
        '<PY:END_YEAR>': end_year,
        '<PY:WINDOW_SIZE>': window_size,
        '<PY:WINDOW_SIZE_GP>': window_size_gp,
        '<PY:CENTER_MONTH>': center_month,
        '<PY:N_PARPOOL>': '8',
        '<PY:OCEAN_BASIN>': ob,
        '<PY:OB_MESHGRID>': ob_mesh,
    })
