#!/usr/bin/env python

import csands_folds_reader,math

#################################

def smoothSpeed(speeds):
	count  = 0
	total = 0
	for speed in speeds:
		if speed:
			count += 1
			total += speed
	if count > (len(speeds)-1)/2:
		return total / count
	return None

def smoothAll(speeds, window):
	output = []
	for i in xrange(len(speeds)):
		if speeds[i]:
			j = max(i - (window - 1) / 2, 0)
			k = min(i + (window + 1) / 2, len(speeds))
			output.append(smoothSpeed(speeds[j:k]))
		else:
			output.append(None)
	return output

#########################################

def get_optimality_thresholds(organisms,scores,window_size,low_thresh,high_thresh,reading_frame_file):

	speed_obs = dict()
	for org in organisms:
		speed_obs[org] = []

	data = csands_folds_reader.get_data(reading_frame_file)
	for datum in data:
		if datum.organism in organisms:
			speed_profile = []
			for i in xrange(len(datum.rna_sequence)/3):
				codon = datum.rna_sequence[i*3:(i+1)*3]
				#try statement to deal with non-standard nucleotides
				try:
					speed_profile.append(scores[datum.organism][codon])
				except:
					continue
			speed_obs[datum.organism].extend(smoothAll(speed_profile,window_size))

	speed_low = dict()
	speed_high = dict()
	for org in speed_obs.keys():
		speed_obs[org].sort()
		speed_low[org] = speed_obs[org][int(math.floor(low_thresh*len(speed_obs[org])))]
		speed_high[org] = speed_obs[org][int(math.ceil(high_thresh*len(speed_obs[org])))]

	return speed_low,speed_high

############################################

def score_rna(aligned_rna,scores):
	speed_profile = []
	for i in xrange(len(rna_seq)/3):
		codon = rna_seq[i*3:(i+1)*3]
		speed_profile.append(scores[codon])
	return speed_profile

def score_rna_multiple(aligned_rnas,organisms,scores):
	speed_profiles = []
	for rna,org in zip(aligned_rnas,organism):
		speed_profiles.append(score_rna(rna,org,scores[org]))
	return speed_profiles
	

