#!/usr/bin/env python
from copy import deepcopy
import wave
import struct
import sys
import numpy as np
from math import sqrt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab
import matplotlib.pyplot as plt
import os
from subprocess import check_call
from tempfile import mktemp
from scikits.audiolab import wavread, play
from scipy.signal import remez, lfilter
from pylab import *
import json

final_val = []
final_val_slice = []

VIOLET = "#ee82ee"
BLUE = "#0000ff"
GREEN = "#008000"
YELLOW = "#ffff00"
ORANGE = "#ffa500"
RED = "#ff0000"

color = [VIOLET, BLUE, GREEN, YELLOW, ORANGE, RED]

def mp3_to_wave(mp3_filename):
	# convert mp3, read wav
	wname = mp3_filename[:-3] + "wav"
	check_call(['avconv', '-i', mp3_filename, wname])
	return wname

def get_color_hash(max_pos):
	count = 0
	inc = 2
	while max_pos > inc:
		inc += 2
		count += 1
	return color[count]

def process_final_val(number_of_peers):
	max_amp = 0
	max_amp = np.amax(final_val)
	#print final_val
	#raw_input()
	print type(final_val)
	print final_val
	final_val_norm = []
	#raw_input()
	for ts in final_val :
		tmp = []
		for i in ts : 
			tmp.append(i/max_amp)
		final_val_norm.append(deepcopy(tmp))
	print final_val_norm

	#final_val_norm = final_val/max_amp
	#print final_val_norm
	#raw_input()
	#max_val = np.amax(final_val_norm)
	#response_array = [[] for i in range(number_of_peers)]
	response_array = []
	for i in final_val_norm:
		for time_slice in i:
			print time_slice
#			raw_input()
	#	print time_slice
	#	raw_input()
			print time_slice
			if time_slice == [] :
				continue
			max_amp_val_ts = max(time_slice)
	#	print max_amp_val_ts
	#	raw_input()
			max_pos = [i for i, j in enumerate(time_slice) if j == max_amp_val_ts]
	#	print max_pos
	#	raw_input()
		#for peer in range(0, len(time_slice)) :
			resp = {}
			resp["Intensity"] = max_amp_val_ts
			resp["Color"] = min(max_pos)
	#	print resp 
	#	raw_input()
			response_array.append(resp)
	  	#response_array[peer].append(resp)
	return response_array


def compute(filename, number_of_peers):
	# Open the wave file and get info
	wave_file_name = mp3_to_wave(filename)
	wave_file = wave.open(wave_file_name, 'r')
	data_size = wave_file.getnframes()
	sample_rate = wave_file.getframerate()
	#print "Sample rate: %s" % sample_rate
	sample_width = wave_file.getsampwidth()
	duration = data_size / float(sample_rate)

	# Read in sample data
	sound_data = wave_file.readframes(data_size)

	#print type(sound_data)
	#print len(sound_data)

	# Close the file, as we don't need it any more
	wave_file.close()

	# Unpack the binary data into an array
	unpack_fmt = '%dh' % (data_size)
	sound_data = struct.unpack(2*unpack_fmt, sound_data)

	# Process many samples
	fouriers_per_second = 8 # Frames per second
	fourier_spread = 1.0/fouriers_per_second
	fourier_width = fourier_spread
	fourier_width_index = fourier_width * float(sample_rate)

	if len(sys.argv) < 3:
		length_to_process = int(duration)-1
	else:
		length_to_process = float(sys.argv[2])

	#print "Fourier width: %s" % str(fourier_width)

	total_transforms = int(round(length_to_process * fouriers_per_second))
	fourier_spacing = round(fourier_spread * float(sample_rate))

	print "Duration: %s" % duration
	print "For Fourier width of "+str(fourier_width)+" need "+str(fourier_width_index)+" samples each FFT"
	print "Doing "+str(fouriers_per_second)+" Fouriers per second"
	print "Total " + str(total_transforms * fourier_spread)
	print "Spacing: "+str(fourier_spacing)
	print "Total transforms "+str(total_transforms)

	lastpoint=int(round(length_to_process*float(sample_rate)+fourier_width_index))-1

	sample_size = fourier_width_index
	freq = sample_rate / sample_size * np.arange(sample_size)

	#x_axis = range(0, 12)

	def getBandWidth():
		return (2.0/sample_size) * (sample_rate / 2.0)

	def freqToIndex(f):
		# If f (frequency is lower than the bandwidth of spectrum[0]
		if f < getBandWidth()/2:
			return 0
		if f > (sample_rate / 2) - (getBandWidth() / 2):
			return sample_size -1
		fraction = float(f) / float(sample_rate)
		index = round(sample_size * fraction)
		return index

	fft_averages = []	

	def average_fft_bands(fft_array):
		num_bands = number_of_peers # The number of frequency bands (12 = 1 octave)
		#final_val += fft_averages
		#print "-----------------------------------------------------------------------------"
		#print "before"
		#print final_val_slice 
		#raw_input()
		#print fft_averages
		#raw_input()
		#print "after"
		final_val.append(deepcopy(fft_averages))
		#print final_val_slice
		#raw_input()
		del fft_averages[:]
		for band in range(0, num_bands):
			avg = 0.0

			if band == 0:
				lowFreq = int(0)
			else:
				lowFreq = int(int(sample_rate / 2) / float(2 ** (num_bands - band)))
			hiFreq = int((sample_rate / 2) / float(2 ** ((num_bands-1) - band)))
			lowBound = int(freqToIndex(lowFreq))
			hiBound = int(freqToIndex(hiFreq))
			for j in range(lowBound, hiBound):
				avg += fft_array[j]

			avg /= (hiBound - lowBound + 1)
			fft_averages.append(avg)
		#print fft_averages
		#raw_input()
		return fft_averages

	for offset in range(0, total_transforms):
		start = int(offset * sample_size)
		end = int((offset * sample_size) + sample_size -1)

#		print "Processing sample %i of %i (%d seconds)" % (offset + 1, total_transforms, end/float(sample_rate))
		sample_range = sound_data[start:end]
		## FFT the data
		fft_data = abs(np.fft.fft(sample_range))
		# Normalise the data a second time, to make numbers sensible
		fft_data *= ((2**.5)/sample_size)
		#plt.ylim(0, 1000)
		average_fft_bands(fft_data)
		#print averages
		#raw_input()
		#final_val.append(averages)
		#averages = []
		#print final_val
		#raw_input()
		#y_axis = fft_averages
		"""Stuff for bar graph"""
		#width = 0.35
		#p1 = plt.bar(x_axis, y_axis, width, color='r')
		"""End bar graph stuff"""
		#filename = str('frame_%05d' % offset) + '.png'
		#plt.savefig(filename, dpi=100)
		#plt.close()
		"""End bar graph stuff"""
		#filename = str('frame_%05d' % offset) + '.png'
		#plt.savefig(filename, dpi=100)
	#print final_val
	return process_final_val(number_of_peers)

if __name__ == '__main__':
	compute(sys.argv[1], 12)
