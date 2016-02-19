#!/usr/bin/env python

import sys 
import yaml
import numpy
import pyfits
import scipy
from scipy.stats import chi2

#import ROOT 

__version__ = 1.0

def get_offset_time_array(fitsfile,extname='EVENTS',time0=None):
	print "fitsfile:%s" % fitsfile
	hdu = pyfits.open(fitsfile)
	time_array = hdu[extname].data['TIME']
	if time0 is None:
		time0 = time_array[0]
	offset_time_array = time_array - time0
	return offset_time_array

def convert_to_phase_array(time_array, period, pdot):
	# ----- Reference of the calcuration -----
	# When constant P and Pdot are given, phases when each X-ray photons
	# comes are calcurated by integrtion of 1/(P + t*Pdot) from 0 to t.
	# Here t is a difference in time of the arrival time and the epoch.
	# This gives,
	#     1/(Pdot){ ln(1 + Pdot/P *t) }
	# ----------------------------------------
	if pdot == 0.0:
		phase_array = (time_array / period) % 1.0
	elif pdot != 0.0:
		if time_array[0] < 0.0:
			print "If pdot!=0, please input epoch t0 smaller than the time of all the event"
			quit() 
		phase_array = numpy.log(1.0 + pdot * time_array / period)/pdot % 1.0

	"""
	hist = ROOT.TH1D("hist","hist",100,0.0,1.0)		
	for phase in phase_array:
		hist.Fill(phase)
	can = ROOT.TCanvas('can','can',800,600)
	hist.Draw('hist')
	can.Print('test.pdf')
	"""	

	return phase_array

def get_folding_statistics(phase_array, nbin, flagHistgramOut=False):
 		"""
 		Standard Folding Analysis. Return Statistics:
 		Chi2-value, reduced-chi2, and survival statistics.
 		"""
 		(count, bins) = numpy.histogram(phase_array, bins=nbin, range=(0.0,1.0))
 		x = 0.5*(bins[1:]+bins[:-1])
 		y = count 
 		x_err = 0.5 * (bins[1:]-bins[:-1])
 		y_err = numpy.sqrt(count)
 		const = scipy.polyfit(x,y,0)[0]
 		chi2 = numpy.sum((y-const)**2/const)
 		rchi2 = chi2/(nbin-1)
 		sf = scipy.stats.chi2.sf(chi2, nbin-1)
 		if flagHistgramOut == True:
 			return chi2, rchi2, sf, [x, y, x_err, y_err]
 		else:
	 		return chi2, rchi2, sf

def get_Zn_statistics(phase_array, max_n=10):
	"""
	Zn-test statistics. See Section 3 Brazier (1994)
	http://adsabs.harvard.edu/abs/1994MNRAS.268..709B
	z2_array ==> Z_n^2 values 
	"""
	number_of_events = len(phase_array)
	cos_square_array = numpy.array([numpy.sum(numpy.cos(j * 2*numpy.pi*phase_array))**2 for j in range(1,max_n+1)])
	sin_square_array = numpy.array([numpy.sum(numpy.sin(j * 2*numpy.pi*phase_array))**2 for j in range(1,max_n+1)])
	z2_array = numpy.array( [ (numpy.sum(cos_square_array[0:j]) + numpy.sum(sin_square_array[0:j]))*2/number_of_events for j in range(1,max_n+1)]  )
	sf_array = numpy.array( [ scipy.stats.chi2.sf(z2_array[j], 2*(j+1)) for j in range(max_n)])
	statistics_array = [ [z2_array[j], sf_array[j]] for j in range(max_n)]
	return statistics_array

def run_pulse_search(yamlfile):
	f = open(yamlfile)
	config = yaml.load(f)
	f.close()

	f = open('%s/%s' % (config['outdir_path'],config['outname']),'w')

	offset_time_array = get_offset_time_array(
		config['input_fitsfile'],
		time0=config['epoch'])

	for i in range(config['search_number']):
		period  = config['period']
		period += (float(i)-0.5*float(config['search_number']))*float(config['search_period_resoltion'])
		phase_array = convert_to_phase_array(offset_time_array, 
			period=period, pdot=config['pdot'])
		out  = '%.17f  ' % period
		chi2, rchi2, sf = get_folding_statistics(phase_array, nbin=100, flagHistgramOut=False)
		out += '%.3e %.3e %.3e  ' % (chi2, rchi2, sf)
		statistics_array = get_Zn_statistics(phase_array, max_n=3)
		for j in statistics_array:
			out += '%d %.3e %.3e  ' % (statistics_array.index(j), j[0],j[1])
		out += '\n'
		f.write(out)
	f.close()

if __name__ == "__main__":

	if len(sys.argv) != 2:
		sys.stderr.write('%s input.yaml \n' % sys.argv[0])
		quit()
	yamlfile    = sys.argv[1]
	run_pulse_search(yamlfile)



