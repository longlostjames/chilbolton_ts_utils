#!/usr/bin/env python
# coding: utf-8

import sys, getopt

import wivern_chilbolton_utils as wivern
import os, getpass, glob

user = getpass.getuser()

import netCDF4 as nc4
import pyart
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import cmocean
import glob
import numpy as np;
import numpy.ma as ma
import cftime
import fnmatch

from datetime import datetime

matplotlib.use('Agg');

from mpl_toolkits.axes_grid1 import make_axes_locatable


def load_wivern2_l1(file):

    """This routine loads a Level 1 WIVERN-2 I/Q time-series NetCDF file.
    It then creates pulse-pair estimates of the Doppler velocity and radar
    equivalent reflectivity factor, annd stores the results in a Python
    dictionary for further use.  The velocity data may be folded and no attempt
    at unfolding is made.

    :param file: Full path to Level 1 file
    :type file: str
    """

    DS = nc4.Dataset(file);

    dtime = cftime.num2pydate(DS['time'],DS['time'].units);

    I = ma.masked_where(DS['qc_flag'][:,:,:]>2,DS['I'][:,:,:]);
    Q = ma.masked_where(DS['qc_flag'][:,:,:]>2,DS['Q'][:,:,:]);

    range_km = DS['range'][:]/1000.;
    altitude_m = DS['altitude'][:];

    nray,npulse,ngate = I.shape;

    # Pulse-pair reflectivity (no distinction made between polarizations)
    # -------------------------------------------------------------------
    Zpp = np.empty([nray,ngate],dtype='float');
    Zpp = np.mean(I**2+Q**2,axis=1);
    Zpp = 10*np.log10(Zpp)+20*np.log10((range_km[None,:]));

    I0  = I[:,2::2,:];
    Q0  = Q[:,2::2,:];
    Im1 = I[:,1:-1:2,:];
    Qm1 = Q[:,1:-1:2,:];
    Im2 = I[:,:-2:2,:];
    Qm2 = Q[:,:-2:2,:];
    Ip1 = I[:,3::2,:];
    Qp1 = Q[:,3::2,:];

    real_vh_arr = I0*Im2+Q0*Qm2;
    imag_vh_arr = Q0*Im2-I0*Qm2;
    real_vv_arr = Ip1*Im1+Qp1*Qm1;
    imag_vv_arr = Qp1*Im1-Ip1*Qm1;

    real_vh = np.sum(real_vh_arr,axis=1);
    imag_vh = np.sum(imag_vh_arr,axis=1);
    real_vv = np.sum(real_vv_arr,axis=1);
    imag_vv = np.sum(imag_vv_arr,axis=1);

    prf  = DS['prf'][:];
    freq = DS['frequency'][:]*1e9;  # convert Ghz to Hz
    c = 299792458; # m/s
    wavelength = c/freq;
    vf = wavelength*prf/4;  # folding velocity

    # Determine pulse-pair velocity
    # -----------------------------
    velpp=ma.arctan2(imag_vh+imag_vv,real_vh+real_vv)*vf/np.pi;

    data = dict();
    data['I'] = I;
    data['Q'] = Q;
    data['Zpp'] = Zpp;
    data['velpp'] = velpp;
    data['range'] = range_km;
    data['altitude'] = altitude_m;
    data['datetime'] = dtime;
    data['folding_velocity'] = vf;

    DS.close();

    return data


def create_wivern2_l1_plot(files,figpath,plotfile,dt_min,dt_max,time_alignment='post'):

    """This routine accepts a list of Level 1 WIVERN-2 I/Q time-series files,
    and generates a composite quicklook plot of radar equivalent reflectivity
    factor and Doppler velocity.

    :param files: List of files to be processed
    :type infile: str

    :param figpath: path where quicklook file is to be written
    :type figpath: str

    :param plotfile: name of quicklook file
    :type plotfile: str

    :param dt_min: lower limit for time axis
    :type dt_min: python datetime

    :param dt_max: upper limit for time axis
    :type dt_max: python datetime

    :param time_alignment: 'post' (default) specifies that times are aligned
        with the end of each ray; 'pre' specifies that times are aligned with the
        start of each ray (appropriate for 35 GHz data);
    :type time_alignment: str
    """

    DS = nc4.Dataset(files[0]);

    radar = DS.source;

    DS.close()

    infile = files[0];
    nfiles = len(files);

    ifile = 1;
    print('{}/{} {}'.format(ifile,nfiles,infile));


    DS0 = load_wivern2_l1(infile);

    vf = DS0['folding_velocity'];

    fig, ax = plt.subplots(2,1,figsize=(8,8),constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.1,wspace=0.2)


    fig.patch.set_facecolor('white')

    myFmt = mdates.DateFormatter('%H:%M')

    ax[0].xaxis.set_major_formatter(myFmt);
    ax[1].xaxis.set_major_formatter(myFmt);

    hmax = 12;

    ax[0].set_xlim(dt_min,dt_max);
    ax[0].set_ylim(0,hmax);
    ax[1].set_xlim(dt_min,dt_max);
    ax[1].set_ylim(0,hmax);

    ax[0].set_xlabel('Time (UTC)');
    ax[1].set_xlabel('Time (UTC)');

#    ax[0].set_ylabel('Height above radar (km)');
#    ax[1].set_ylabel('Height above radar (km)');

    ax[0].set_ylabel('Altitude (km)');
    ax[1].set_ylabel('Altitude (km)');

    gate_width = DS0['range'][1]-DS0['range'][0];
    gate_starts = DS0['range'][:] - gate_width/2.0;
    gate_starts = np.append(gate_starts,DS0['range'][-1]+gate_width/2.0);
    gate_starts = gate_starts + DS0['altitude']/1000.;

    timestep = DS0['datetime'][1]-DS0['datetime'][0];
    dt = DS0['datetime'][:];
    dt_pre = DS0['datetime'][0]-timestep;
    dt_post = DS0['datetime'][-1]+timestep;

    if time_alignment=='post':
        dt = np.insert(dt,[0],dt_pre);
    elif time_alignment=='pre':
        dt = np.append(dt,dt_post);

    vel_cmap = plt.cm.get_cmap('RdBu')
    vel_cmap = vel_cmap.reversed()

    h1 = ax[0].pcolor(dt,gate_starts,DS0['Zpp'][:,:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
    cb1 = plt.colorbar(h1,ax=ax[0],orientation='horizontal',shrink=0.8);
    cb1.ax.set_xlabel("Reflectivity Factor (dBZ)");

    titlestr = radar;
    ax[0].set_title(titlestr,loc='left')
    ax[0].set_title('{}'.format(dt_min.date()),loc='right')

    ax[1].xaxis.set_major_formatter(myFmt);
    h2 = ax[1].pcolor(dt,gate_starts,DS0['velpp'][:,:].transpose(),vmin=-vf,vmax=vf,cmap=vel_cmap);
    cb2 = plt.colorbar(h2,ax=ax[1],orientation='horizontal',shrink=0.8);
    cb2.ax.set_xlabel("Doppler velocity (m/s)");

    titlestr = radar;
    ax[1].set_title(titlestr,loc='left')
    ax[1].set_title('{}'.format(dt_min.date()),loc='right')

    for file in files[1:]:

        ifile = ifile+1;

        print('{}/{} {}'.format(ifile,nfiles,file));

        DS = nc4.Dataset(file);
        DS_start_last = cftime.num2pydate(DS['time'][-1],DS['time'].units);
        DS_end_first = cftime.num2pydate(DS['time'][0],DS['time'].units);
        DS_start_first = cftime.num2pydate(DS['time'][0],DS['time'].units);
        DS_end_last = cftime.num2pydate(DS['time'][-1],DS['time'].units);
        DS.close();

        if ((DS_start_last >= dt_min) & (DS_start_first <= dt_max)) | ((DS_end_last >= dt_min) & (DS_end_first <= dt_max)):
            DS0 = load_wivern2_l1(file);

            gate_width = DS0['range'][1]-DS0['range'][0];
            gate_starts = DS0['range'][:] - gate_width/2.0;
            gate_starts = np.append(gate_starts,DS0['range'][-1]+gate_width/2.0);
            gate_starts = gate_starts + DS0['altitude']/1000.;

            timestep = DS0['datetime'][1]-DS0['datetime'][0];
            dt = DS0['datetime'][:];
            dt_pre = DS0['datetime'][0]-timestep;
            dt_post = DS0['datetime'][-1]+timestep;

            if time_alignment=='post':
                dt = np.insert(dt,[0],dt_pre);
            elif time_alignment=='pre':
                dt = np.append(dt,dt_post);

            ax[0].pcolor(dt,gate_starts,DS0['Zpp'][:,:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
            ax[1].pcolor(dt,gate_starts,DS0['velpp'][:,:].transpose(),vmin=-vf,vmax=vf,cmap=vel_cmap);


            ax[0].grid(True);
            ax[1].grid(True);

#    plt.tight_layout(w_pad=2 / 72, h_pad=2 / 72);
#    plt.subplots_adjust(hspace=0.4,wspace=0);


    plt.savefig(os.path.join(figpath,plotfile));

    plt.show();
    plt.close();

basepath = '/gws/nopw/j04/ncas_obs/cao/processing/';

camra_l1path = '/gws/nopw/j04/ncas_obs/cao/processing/ncas-radar-camra-1/campaign/wivern2/L1';
copernicus_l1path = '/gws/nopw/j04/ncas_obs/cao/processing/ncas-radar-ka-band-1/campaign/wivern2/L1';
galileo_l1path = '/gws/nopw/j04/ncas_obs/cao/processing/ncas-radar-w-band-1/campaign/wivern2/L1';
figpath = '/home/users/cjwalden/jupyter/wivern2'

def main(argv):

   outpath = figpath;

   try:
      opts, args = getopt.getopt(argv,"hr:d:i:o:m:",["radar=","date=","inpath=","outpath="])
   except getopt.GetoptError:
      print ('wivern2_quicklooks.py -r <radar> -d <date> -i <input_path> -o <output_path>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('wivern2_quicklooks.py -r <radar> -d <date> -i <input_path> -o <output_path>')
         sys.exit()
      elif opt in ("-r", "--radar"):
         radar = arg
      elif opt in ("-d", "--date"):
         datestr = arg
      elif opt in ("-i", "--inpath"):
         inpath = arg
      elif opt in ("-o", "--outpath"):
         outpath = arg
      elif opt in ("-m", "--minmax"):
         duration = arg
#   print ('Radar is ', radar);
#   print ('Input path is ', inpath);
#   print ('Output path is ', outpath)

   inpath3 = os.path.join(basepath,'ncas-radar-camra-1','campaign','wivern2','L1',datestr);
   inpath35 = os.path.join(basepath,'ncas-radar-ka-band-1','campaign','wivern2','L1',datestr);
   inpath94 = os.path.join(basepath,'ncas-radar-w-band-1','campaign','wivern2','L1',datestr);

   if radar in ("camra", "CAMRa"):
       ncas_radar = "ncas-radar-camra-1"
       print("Making CAMRa quicklooks for {}".format(datestr));
       make_quicklooks_single(ncas_radar,datestr,basepath,outpath);
   elif radar in ("copernicus", "Copernicus"):
       ncas_radar = "ncas-radar-ka-band-1"
       print("Making Copernicus quicklooks for {}".format(datestr));
       make_quicklooks_single(ncas_radar,datestr,basepath,outpath,time_alignment='pre');
   elif radar in ("galileo", "Galileo"):
       ncas_radar = "ncas-radar-w-band-1"
       print("Making Galileo quicklooks for {}".format(datestr));
       make_quicklooks_single(ncas_radar,datestr,basepath,outpath);
   elif radar in ("multi", "Multi"):
       make_quicklooks_multi(datestr,basepath,outpath,duration);


def make_quicklooks_single(ncas_radar,datestr,basepath,outpath,time_alignment='post'):

    pattern = '*_fix-ts_l1_v1.0.nc'

    inpath = os.path.join(basepath,ncas_radar,'campaign','wivern2','L1',datestr);

    l1dirs  = [];

    for root,dirs,files in os.walk(inpath):
        l1dirs += [os.path.join(root,d) for d in dirs];

    print(l1dirs);

    if l1dirs==[]:
        l1dirs  = inpath;
        l1files  = [os.path.join(inpath,f) for f in fnmatch.filter(get_files(inpath), pattern)];

        l1files.sort();

        dt_start = datetime(3000,1,1);
        dt_end = datetime(1900,1,1);

        if len(l1files)>0:
            DSstart  = nc4.Dataset(l1files[0]);
            dt_start = cftime.num2pydate(DSstart['time'][0],DSstart['time'].units);
            DSend    = nc4.Dataset(l1files[-1]);
            dt_end   = cftime.num2pydate(DSend['time'][-1],DSend['time'].units);
            DSstart.close();
            DSend.close();

        print("Start = {}".format(dt_start));
        print("End = {}".format(dt_end));

        if len(l1files)>0:
            plotfile  = os.path.join('{}_esa-wivern2_{}.png'.format(ncas_radar,datestr));
            create_wivern2_l1_plot(l1files,outpath,plotfile,dt_start,dt_end,time_alignment);
    else:
        for d in l1dirs:
            print(d);
            event = os.path.split(d)[-1];
            print(event);

            dir = os.path.join(inpath,event);
            l1files  = [os.path.join(inpath,f) for f in fnmatch.filter(get_files(dir), pattern)];

            l1files.sort();

            dt_start = datetime(3000,1,1);
            dt_end = datetime(1900,1,1);

            if len(l1files)>0:
                DSstart  = nc4.Dataset(l1files[0]);
                dt_start = cftime.num2pydate(DSstart['time'][0],DSstart['time'].units);
                DSend    = nc4.Dataset(l1files[-1]);
                dt_end   = cftime.num2pydate(DSend['time'][-1],DSend['time'].units);
                DSstart.close();
                DSend.close();

            print("Start = {}".format(dt_start));
            print("End = {}".format(dt_end));


            if len(l1files)>0:
                plotfile  = os.path.join('{}_esa-wivern2_{}_{}.png'.format(ncas_radar,datestr,event));
                create_wivern2_l1_plot(l1files,outpath,plotfile,dt_start,dt_end,time_alignment);

get_files = lambda path: (os.path.join(root, file) for root, dirs, files in os.walk(path) for file in files)

def make_quicklooks_multi(datestr,basepath,outpath,duration='max'):

    """This routine creates a set of time-height quicklooks for WIVERN-2 campaign
    data from a given date.  It assumes a directory tree in which any separate
    events are nested in subdirectories (e.g. EventA, EventB).  The limits of
    the time axis are set to be common for all radars.

    :param datestr: String reresenting date in the format YYYYmmdd
    :type datestr: str

    :param basepath: path where individual radar data directory trees are located
    :type basepath: str

    :param outpath: destination path for quicklooks
    :type outpath: str

    :param duration: 'max' (default) sets the time limits to the maximum across
        all radars, 'min' sets the time limits to the common overlap period.
    :type duration: str
    """

    inpath3  = os.path.join(basepath,'ncas-radar-camra-1','campaign','wivern2','L1',datestr);
    inpath35 = os.path.join(basepath,'ncas-radar-ka-band-1','campaign','wivern2','L1',datestr);
    inpath94 = os.path.join(basepath,'ncas-radar-w-band-1','campaign','wivern2','L1',datestr);

    pattern = '*_fix-ts_l1_v1.0.nc'

    l1dirs3  = [];
    l1dirs35 = [];
    l1dirs94 = [];

    for root,dirs,files in os.walk(inpath3):
        l1dirs3 += [os.path.join(root,d) for d in dirs];

    for root,dirs,files in os.walk(inpath35):
        l1dirs35 += [os.path.join(root,d) for d in dirs];

    for root,dirs,files in os.walk(inpath94):
        l1dirs94 += [os.path.join(root,d) for d in dirs];

    print(l1dirs3);
    print(l1dirs35);
    print(l1dirs94);

    if l1dirs3==[]:
        l1dirs3  = inpath3;
        l1dirs35 = inpath35;
        l1dirs94 = inpath94;
        l1files3  = [os.path.join(inpath3,f) for f in fnmatch.filter(get_files(inpath3), pattern)];
        l1files35 = [os.path.join(inpath35,f) for f in fnmatch.filter(get_files(inpath35), pattern)];
        l1files94 = [os.path.join(inpath94,f) for f in fnmatch.filter(get_files(inpath94), pattern)];

        l1files3.sort();
        l1files35.sort();
        l1files94.sort();

        if duration=='max':
            dt_start = datetime(3000,1,1);
            dt_end = datetime(1900,1,1);
        else:
            dt_start = datetime(1900,1,1);
            dt_end = datetime(3000,1,1);

        if len(l1files3)>0:
            DS3start  = nc4.Dataset(l1files3[0]);

            timestep = DS3start['time'][1]-DS3start['time'][0];

            if duration=='max':
                dt_start   = min(dt_start,cftime.num2pydate(DS3start['time'][0]-timestep,DS3start['time'].units));
            else:
                dt_start   = max(dt_start,cftime.num2pydate(DS3start['time'][0]-timestep,DS3start['time'].units));
            DS3end    = nc4.Dataset(l1files3[-1]);
            if duration=='max':
                dt_end     = max(dt_end,cftime.num2pydate(DS3end['time'][-1],DS3end['time'].units));
            else:
                dt_end     = min(dt_end,cftime.num2pydate(DS3end['time'][-1],DS3end['time'].units));
            DS3start.close();
            DS3end.close();

        if len(l1files35)>0:
            DS35start  = nc4.Dataset(l1files35[0]);

            if duration=='max':
                dt_start   = min(dt_start,cftime.num2pydate(DS35start['time'][0],DS35start['time'].units));
            else:
                dt_start   = max(dt_start,cftime.num2pydate(DS35start['time'][0],DS35start['time'].units));
            DS35end    = nc4.Dataset(l1files35[-1]);
            timestep = DS35end['time'][-1]-DS35end['time'][-2];

            if duration=='max':
                dt_end     = max(dt_end,cftime.num2pydate(DS35end['time'][-1]+timestep,DS35end['time'].units));
            else:
                dt_end     = min(dt_end,cftime.num2pydate(DS35end['time'][-1]+timestep,DS35end['time'].units));
            DS35start.close();
            DS35end.close();

        if len(l1files94)>0:
            DS94start  = nc4.Dataset(l1files94[0]);
            timestep = DS94start['time'][1]-DS94start['time'][0];

            if duration=='max':
                dt_start   = min(dt_start,cftime.num2pydate(DS94start['time'][0]-timestep,DS94start['time'].units));
            else:
                dt_start   = max(dt_start,cftime.num2pydate(DS94start['time'][0]-timestep,DS94start['time'].units));
            DS94end    = nc4.Dataset(l1files94[-1]);
            if duration=='max':
                dt_end     = max(dt_end,cftime.num2pydate(DS94end['time'][-1],DS94end['time'].units));
            else:
                dt_end     = min(dt_end,cftime.num2pydate(DS94end['time'][-1],DS94end['time'].units));
            DS94start.close();
            DS94end.close();


        print("Start = {}".format(dt_start));
        print("End = {}".format(dt_end));


        if len(l1files3)>0:
            plotfile3  = os.path.join('ncas-radar-camra-1_esa-wivern2_{}.png'.format(datestr));
            create_wivern2_l1_plot(l1files3,outpath,plotfile3,dt_start,dt_end);
        if len(l1files35)>0:
            plotfile35 = os.path.join('ncas-radar-ka-band-1_esa-wivern2_{}.png'.format(datestr));
            create_wivern2_l1_plot(l1files35,outpath,plotfile35,dt_start,dt_end,time_alignment='pre');
        if len(l1files94)>0:
            plotfile94 = os.path.join('ncas-radar-w-band-1_esa-wivern2_{}.png'.format(datestr));
            create_wivern2_l1_plot(l1files94,outpath,plotfile94,dt_start,dt_end);
    else:
        for d in l1dirs3:
            print(d);
            event = os.path.split(d)[-1];
            print(event);
            d3  = d;
            d35 = os.path.join(inpath35,event);
            d94 = os.path.join(inpath94,event);
            l1files3  = [os.path.join(inpath3,f) for f in fnmatch.filter(get_files(d3), pattern)];
            l1files35 = [os.path.join(inpath35,f) for f in fnmatch.filter(get_files(d35), pattern)];
            l1files94 = [os.path.join(inpath94,f) for f in fnmatch.filter(get_files(d94), pattern)];

            l1files3.sort();
            l1files35.sort();
            l1files94.sort();

            if duration=='max':
                dt_start = datetime(3000,1,1);
                dt_end = datetime(1900,1,1);
            else:
                dt_start = datetime(1900,1,1);
                dt_end = datetime(3000,1,1);

            if len(l1files3)>0:
                DS3start  = nc4.Dataset(l1files3[0]);
                timestep = DS3start['time'][1]-DS3start['time'][0];

                if duration=='max':
                    dt_start   = min(dt_start,cftime.num2pydate(DS3start['time'][0]-timestep,DS3start['time'].units));
                else:
                    dt_start   = max(dt_start,cftime.num2pydate(DS3start['time'][0]-timestep,DS3start['time'].units));
                DS3end    = nc4.Dataset(l1files3[-1]);
                if duration=='max':
                    dt_end     = max(dt_end,cftime.num2pydate(DS3end['time'][-1],DS3end['time'].units));
                else:
                    dt_end     = min(dt_end,cftime.num2pydate(DS3end['time'][-1],DS3end['time'].units));
                DS3start.close();
                DS3end.close();

            if len(l1files35)>0:
                DS35start  = nc4.Dataset(l1files35[0]);
                if duration=='max':
                    dt_start   = min(dt_start,cftime.num2pydate(DS35start['time'][0],DS35start['time'].units));
                else:
                    dt_start   = max(dt_start,cftime.num2pydate(DS35start['time'][0],DS35start['time'].units));
                DS35end    = nc4.Dataset(l1files35[-1]);
                timestep = DS35end['time'][-1]-DS35end['time'][-2];
                if duration=='max':
                    dt_end     = max(dt_end,cftime.num2pydate(DS35end['time'][-1]+timestep,DS35end['time'].units));
                else:
                    dt_end     = min(dt_end,cftime.num2pydate(DS35end['time'][-1]+timestep,DS35end['time'].units));
                DS35start.close();
                DS35end.close();

            if len(l1files94)>0:
                DS94start  = nc4.Dataset(l1files94[0]);
                timestep = DS94start['time'][1]-DS94start['time'][0];
                if duration=='max':
                    dt_start   = min(dt_start,cftime.num2pydate(DS94start['time'][0]-timestep,DS94start['time'].units));
                else:
                    dt_start   = max(dt_start,cftime.num2pydate(DS94start['time'][0]-timestep,DS94start['time'].units));
                DS94end    = nc4.Dataset(l1files94[-1]);
                if duration=='max':
                    dt_end     = max(dt_end,cftime.num2pydate(DS94end['time'][-1],DS94end['time'].units));
                else:
                    dt_end     = min(dt_end,cftime.num2pydate(DS94end['time'][-1],DS94end['time'].units));
                DS94start.close();
                DS94end.close();

            print("Start = {}".format(dt_start));
            print("End = {}".format(dt_end));


            if len(l1files3)>0:
                plotfile3  = os.path.join('ncas-radar-camra-1_esa-wivern2_{}_{}.png'.format(datestr,event));
                create_wivern2_l1_plot(l1files3,outpath,plotfile3,dt_start,dt_end);
            if len(l1files35)>0:
                plotfile35 = os.path.join('ncas-radar-ka-band-1_esa-wivern2_{}_{}.png'.format(datestr,event));
                create_wivern2_l1_plot(l1files35,outpath,plotfile35,dt_start,dt_end,time_alignment='post');
            if len(l1files94)>0:
                plotfile94 = os.path.join('ncas-radar-w-band-1_esa-wivern2_{}_{}.png'.format(datestr,event));
                create_wivern2_l1_plot(l1files94,outpath,plotfile94,dt_start,dt_end);


if __name__ == "__main__":
   main(sys.argv[1:])
