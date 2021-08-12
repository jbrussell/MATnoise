# %% markdown
# ## miniseed to SAC and remove response
#
# Reads local miniseed data chops, it into day-long files, and converts it to SAC. 
# This version removes instrument response using RESP files
#
# Processing steps:
# - Demean
# - Detrend
# - Trim to 24 hour segment
# - Cosine taper
# - Remove instrument response
# - Anti-alias: low pass (corner = 0.4*sr_new Hz)
# - Remove daily fluctuations: high pass (corner = 1/60/60 Hz)
# - Decimate to 1 Hz
# - Convert to SAC and add station lat, lon, depth info
#
#
# To get a left handed system that matches IRIS, need:
# - CH0 -> BH2
# - CH1 -> BH1
# - CH2 -> BHZ
# - CH3 -> BDH
#
# %% codecell
from obspy import read, read_inventory
from obspy.io.sac import SACTrace
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
# %% codecell
# Setup paths

path2mseed = "Path/to/input/mseed/files" # path to input seed data
path2sac = "Path/to/output/sac/files"  # output sac path
path2sta = "Path/to/stations/file/"  # station file (station ID, lat, lon, elevation)
path2RESP = "Path/to/resp/files/" # instrument response RESP file

sr_new = 1  # new sample rate in Hz
chs = ["CH0", "CH1", "CH2", "CH3"] # channel names

# %% codecell
# Load station file
inventory = pd.read_csv(path2sta, delimiter='\s+', index_col=False, names=['station', 'stla', 'stlo', 'stel'])
inventory = inventory.set_index('station')
# inventory = inventory.loc[station].reset_index()

# %% codecell
# Loop over stations
for ista,sta in enumerate(inventory.station):
    if not os.path.exists(path2sac+sta):
        os.makedirs(path2sac+sta)
    path2mseedsta = path2mseed + sta + '/Mseed/'
    # Loop over channels
    for ch in chs:
        # Get component name (flip BH2 and BH1)
        if ch == "CH0":
            comp = "BH2"
        elif ch == "CH1":
            comp = "BH1"
        elif ch == "CH2":
            comp = "BHZ"
        elif ch == "CH3":
            comp = "BDH"
        pathlist = sorted(Path(path2mseedsta).glob('**/*'+ch+'*.msd'))

        # Loop over day files and extract days
        for ifil, path in enumerate(pathlist):
            if (ifil == 0) or (ifil == len(pathlist)-1):
                continue

            # Load day before, day of, and day after for merging
            st = read(str(pathlist[ifil-1]))
            st += read(str(pathlist[ifil]))
            st += read(str(pathlist[ifil+1]))

            # Check if day already processed
            if len(glob.glob(path2sac+sta+'/'+sta+'.'+str(st[1].stats.starttime.year)+'.'+'%03i'%(st[1].stats.starttime.julday)+'.*.'+comp+'.sac')) != 0:
                print('Skipping... already processed  ' + str(path))
                continue
            print("working on: " + str(path))

            # Check that data segments are neighboring
            sr = st[0].stats.sampling_rate
            dt_strt_gap = st[1].stats.starttime-st[0].stats.endtime
            dt_end_gap = st[2].stats.starttime-st[1].stats.endtime
            if abs(dt_strt_gap*sr) > 2 or abs(dt_end_gap*sr) > 2:
                print('Days are not in sequence! ::: '+str(path))
                print('(Day1 strt) - (Day0 end) = ' + str(st[1].stats.starttime-st[0].stats.endtime) + 'sec')
                print('(Day2 strt) - (Day1 end) = ' + str(st[2].stats.starttime-st[1].stats.endtime) + 'sec')
                continue

            # Get start time as beginning of middle day
            tstart = st[1]
            refdate = st[1].stats.starttime
            tstart = datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0)

            # Merge data (interpolate gaps)
            st.merge(method=1, fill_value='interpolate')

            # Detrend and taper waveform
            st.detrend(type='demean')
            st.detrend(type='linear')

            # make sure trace is 24 hours long
            t1 = UTCDateTime(tstart)
            t2 = t1 + 24*60*60
            st.trim(starttime=t1, endtime=t2, pad=True, nearest_sample=False, fill_value=0)

            # Taper new waveform
            st.taper(type="cosine", max_percentage=0.05)

            # Read RESP file and remove instrument response
            RESPpath = glob.glob(path2RESP+"/"+ch+"/RESP.*."+sta+".."+ch)[0]
            inv = read_inventory(RESPpath)
            sensitivity = inv[0][0][0].response.instrument_sensitivity.value
            resp = inv[0][0][0].response.get_paz()
            poles = resp.poles
            zeros = resp.zeros
            gain = resp.normalization_factor
            units = inv[0][0][0].response.instrument_sensitivity.input_units
            # If units are m/s, add additional zero to convert to m
            if units == "M/S":
                zeros.insert(0,0+0j) # will now be in units of displacement (m)
            # Build dictionary for input to "simulate" function
            paz = {
                'zeros': zeros,
                'poles': poles,
                'gain': gain,
                'sensitivity': sensitivity}
            # Remove instrument response
            st.simulate(paz_remove=paz, paz_simulate=None, pre_filt=[0.001, 0.005, sr/3, sr/2])
#                     st.remove_response(output="DISP", zero_mean=True, taper=True, taper_fraction=0.05, pre_filt=[0.001, 0.005, sr/3, sr/2], water_level=60)

            # Downsample
            st.detrend(type='demean')
            st.detrend(type='linear')
            st.filter('lowpass', freq=0.4*sr_new, zerophase=True)  # anti-alias filter
            st.filter('highpass', freq=1/60/60, zerophase=True)  # Remove daily oscillations
            st.decimate(factor=int(sr/sr_new), no_filter=True)  # downsample
            st.taper(type="cosine", max_percentage=0.05)
            st.detrend(type='demean')
            st.detrend(type='linear')

            # convert to SAC and fill out station/event header info
            sac = SACTrace.from_obspy_trace(st[0])
            sac.stel = inventory.loc[ista].stel
            sac.stla = inventory.loc[ista].stla
            sac.stlo = inventory.loc[ista].stlo
            
            kcmpnm = comp
            sac.kcmpnm = kcmpnm
            yr = str(st[0].stats.starttime.year)
            jday = '%03i'%(st[0].stats.starttime.julday)
            hr = '%02i'%(st[0].stats.starttime.hour)
            mn = '%02i'%(st[0].stats.starttime.minute)
            sec = '%02i'%(st[0].stats.starttime.second)
            sac_out = path2sac + sta + '/' + sta+'.'+yr+'.'+jday+'.'+hr+'.'+mn+'.'+sec+'.'+kcmpnm+'.sac'
            sac.write(sac_out)
