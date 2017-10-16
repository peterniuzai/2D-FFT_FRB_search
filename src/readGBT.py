#! /bin/env python 

import numpy as np
import scipy as sp
from scipy.signal import medfilt
import pyfits
import h5py
import matplotlib.pyplot as plt

def check_fits(data_path):

    hdulist = pyfits.open(data_path, memmap=False)
    print hdulist

    mheader = hdulist[1].header

    for key in mheader.keys():
        print key, '\t', mheader[key]
    #print mheader['STT_IMJD']
    #print mheader['STT_SMJD']
    #print mheader['STT_OFFS']

    mheader = hdulist[0].header

    for key in mheader.keys():
        print key, '\t', mheader[key]
    #print mheader['STT_IMJD']
    #print mheader['STT_SMJD']
    #print mheader['STT_OFFS']

    for k in range(1, 2):
        tbdata = hdulist[k].data

        fieldlabel = []
        for i in range(hdulist[k].header['TFIELDS']): 
            fieldlabel.append(hdulist[k].header['TTYPE%d'%(i+1)])
        print fieldlabel
        #for i in range(len(tbdata)):
        #    print tbdata[i][fieldlabel[k]]
        for i in range(hdulist[k].header['TFIELDS']):
            print hdulist[k].header['TTYPE%d'%(i+1)]
            print sp.unique(tbdata.field(fieldlabel[i])).shape
            print tbdata.field(fieldlabel[i]).shape
            print tbdata.field(fieldlabel[i])
            print

def read_fits(hdulist, nrecords_read=None):

    mheader = hdulist[0].header
    dheader = hdulist[1].header
    freq = hdulist[1].data[0]['DAT_FREQ'].astype(float)
    delta_t = dheader['TBIN']
    nfreq = dheader['NCHAN']
    delta_f = dheader['CHAN_BW']
    freq0 = mheader['OBSFREQ'] - mheader['OBSBW'] / 2. + delta_f / 2
    time0 = mheader["STT_SMJD"] + mheader["STT_OFFS"]

    print "The raw intigration time is %8.5f s"%delta_t

    record_ra  = hdulist[1].data.field('RA_SUB')
    record_dec = hdulist[1].data.field('DEC_SUB')

    nrecords = len(hdulist[1].data)
    if nrecords_read == None:
        nrecords_read = nrecords
    record_ind = np.arange(nrecords_read)

    ntime_record, npol, nfreq, one = hdulist[1].data[0]["DATA"].shape
    print ntime_record, npol, nfreq, one
    out_data = np.empty((nfreq, npol, nrecords_read, ntime_record), dtype=np.float32)

    for ii in range(nrecords_read):
        index = record_ind[ii] 
        #print record_ra[index], record_dec[index]
        record = hdulist[1].data[index]["DATA"]
        scl = hdulist[1].data[index]["DAT_SCL"]
        scl.shape = (1, npol, nfreq, 1)
        offs = hdulist[1].data[index]["DAT_OFFS"]
        offs.shape = (1, npol, nfreq, 1)
        record *= scl
        record += offs
        # Interpret as unsigned int (for Stokes I only).
        record = record.view(dtype=np.uint8)
        # Select stokes I and copy.
        out_data[:,0,ii,:] = np.transpose(record[:,0,:,0])
        # Interpret as signed int (except Stokes I).
        record = record.view(dtype=np.int8)
        out_data[:,1,ii,:] = np.transpose(record[:,1,:,0])
        out_data[:,2,ii,:] = np.transpose(record[:,2,:,0])
        out_data[:,3,ii,:] = np.transpose(record[:,3,:,0])
    out_data.shape = (nfreq, npol, nrecords_read * ntime_record)

    #time = np.arange(data.shape[2]) * delta_t + time0
    time = np.arange(record_ind[0]*ntime_record, (record_ind[-1] + 1)*ntime_record)\
            * delta_t + time0
    subint_time = time0 + hdulist[1].data['OFFS_SUB']
    ra = sample_subint(subint_time, hdulist[1].data['RA_SUB'], time)
    dec = sample_subint(subint_time, hdulist[1].data['DEC_SUB'], time)
    az = sample_subint(subint_time, hdulist[1].data['TEL_AZ'], time)
    el = 90. - sample_subint(subint_time, hdulist[1].data['TEL_ZEN'], time)

    return out_data, time, ra, dec, az, el, freq

def sample_subint(sub_time, sub_var, time):

    diff = np.diff(sub_var) / np.diff(sub_time)
    rate = np.mean(diff)
    start_ind = np.argmin(np.abs(time - sub_time[0]))
    return (time - time[start_ind]) * rate + sub_var[0]

def read_raw(data_path, data_name, nrecords_read=None):

    hdulist = pyfits.open(data_path + data_name + '.fits', 'readonly', memmap=False)
    data, time, ra, dec, az, el, freq = \
            read_fits(hdulist, nrecords_read=nrecords_read)

    timestream_data = TimestreamData(data, time, freq, ra, dec)

    data = remove_cal(timestream_data, p=64, check=False)

    return timestream_data

def write(timestream_data, output_path):

    data_file = h5py.File(output_path + '.hdf5', 'w')
    data_file['data'] = timestream_data.data
    data_file['time'] = timestream_data.time
    data_file['freq'] = timestream_data.freq
    data_file['ra']   = timestream_data.ra
    data_file['dec']  = timestream_data.dec
    #data_file['az']   = az 
    #data_file['el']   = el

    data_file.close()

def load(data_path):

    data_file = h5py.File(data_path + '.hdf5', 'r')

    timestream_data = TimestreamData(
            data_file['data'].value, 
            data_file['time'].value, 
            data_file['freq'].value, 
            data_file['ra'].value, 
            data_file['dec'].value)

    #data_file.close()

    return timestream_data

class TimestreamData(object):


    def __init__(self, data, time, freq, ra, dec):

        self.data = data.view()
        self.time = time.view()
        self.freq = freq.view()
        self.ra   = ra.view()
        self.dec  = dec.view()

def get_phase(data, p=64):

    time_data = np.median(np.mean(data[:, 0, :], axis=0).reshape(-1, p), axis=0)
    noise_level = time_data.max() - time_data.min()
    phase = 0
    while 1:
        if time_data[0] - time_data[-1] > 0.4 * noise_level:
            print phase
            #plt.plot(range(time_data.shape[0]), time_data, '.-')
            #plt.show()
            return phase
        else:
            phase += 1
        time_data = np.roll(time_data, -1, axis=0)

        if phase > p:
            print "error in finding the phase"
            exit()

def remove_cal(data, p=64, check=False):

    phase = get_phase(data.data, p=64)

    index = np.roll(np.arange(64), -1 * phase)
    calon_index = np.sort(index[1:32])
    calof_index = np.sort(index[33:])

    good_index = np.sort(np.concatenate([index[1:32], index[33:]]))

    #print calon_index
    #print calof_index

    shape = data.data.shape
    data.data.shape = shape[:-1] + (shape[-1]/64, 64)
    data.time.shape = (shape[-1]/64, 64)

    #data_calon = data.data[..., calon_index]
    #data_calof = data.data[..., calof_index]

    cal_signal = data.data[..., calon_index] - data.data[..., calof_index]
    cal_signal = np.mean(cal_signal, axis=(-1, -2))

    data.data[..., calon_index] -= cal_signal[..., None, None]
    cal_signal[cal_signal==0] = np.inf
    data.data /= cal_signal[..., None, None]

    data.data = data.data[...,good_index]
    data.time = data.time.reshape(shape[-1]/64, 64)[...,good_index]
    data.ra   = data.ra.reshape(shape[-1]/64, 64)[...,good_index]
    data.dec  = data.dec.reshape(shape[-1]/64, 64)[...,good_index]

    data.time = data.time.flatten()

    data.data = data.data.reshape(
            shape[:-1] + (data.data.shape[-2] * data.data.shape[-1],))

    if check:

        plt.plot(data.time, np.mean(data.data[:, 0, :], axis=0))
        plt.show()


    return data


def check_data(data_path, data_name):

    data = load(data_path + data_name)

    #plt.plot(data.ra, data.dec)
    #plt.show()

    phase = get_phase(data.data, p=64)

    index = np.roll(np.arange(64), -1 * phase)
    calon_index = np.sort(index[1:32])
    calof_index = np.sort(index[33:])

    print calon_index
    print calof_index

    #spec_I = np.mean(data['data'].value[:, 0, 100], axis=0)
    #plt.plot(data['time'].value, spec_I)
    spec_I = data.data[:, 0, :]
    time_I = data.time
    shape = spec_I.shape
    spec_I.shape = shape[:-1] + (shape[-1]/64, 64)
    time_I.shape = (shape[-1]/64, 64)
    spec_I_calon = spec_I[...,calon_index]
    time_I_calon = time_I[...,calon_index].flatten()
    spec_I_calof = spec_I[...,calof_index]
    time_I_calof = time_I[...,calof_index].flatten()

    plt.plot(time_I_calon, np.mean(spec_I_calon, axis=0).flatten())
    plt.plot(time_I_calof, np.mean(spec_I_calof, axis=0).flatten())
    plt.show()

if __name__=="__main__":

    check_fits('/home/ycli/data/gbt/raw/origin/GBT14B_339/78_wigglez1hr_centre_ralongmap_80.fits')

    exit()


    data_path  = '/project/ycli/data/gbt/AGBT14B_339/01_20140718/'
    data_name  = 'guppi_56856_wigglez1hr_centre_0011_0001'

    #output_path = '/project/ycli/data/gbt/converted/'
    output_path = '/project/ycli/data/gbt/cal_removed/'

    check_fits(data_path + data_name + '.fits')

    #data = read_raw(data_path, data_name, 10)
    data = read_raw(data_path, data_name)
    write(data, output_path + data_name)
    #data = load(output_path + data_name)

    #data = remove_cal(data, p=64, check=True)

    #check_data(output_path, data_name)






