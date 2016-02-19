#!/usr/bin/env python

import os 
import sys 
import yaml
import numpy 
import pyfits
from scipy.interpolate import interp1d

__version__ = 1.0

debug = False

TOOLBAR_WIDTH = 40

def randomgen(func,func_range,func_max):
	while(1):
		x_rand, y_rand = numpy.random.random(2)
		x = func_range[0] + (func_range[1]-func_range[0]) * x_rand
		if ( func_max * y_rand <= func(x) ):
			return x
	# http://pacocat.com/?p=596

class PulseProfile():
	def __init__(self, yamlfile):
		sys.stdout.write('----------------------------\n')
		sys.stdout.write('fake_pulse.py version %s\n' % (__version__))
		sys.stdout.write('----------------------------\n')
		if not os.path.exists(yamlfile):
			sys.stderr.write('yamlfile %s does not exists.\n' % yamlfile)
			quit()
		f = open(yamlfile)
		self.config = yaml.load(f)
		f.close()

	def show_parameters(self):
		sys.stdout.write('[show_parameters] %s \n' % self.config)
	
	def read_pulse_profile(self):
		if not os.path.exists(self.config['template_pulse_profile']):
			sys.stderr.write('template_pulse_profile file %s does not exists.\n' % self.config['template_pulse_profile'])
			quit()

		if self.config['flag_template_error']:
			self.phase, self.template_profile, self.template_profile_error = numpy.loadtxt(self.config['template_pulse_profile'],unpack=True)		
		else:
			self.phase, self.template_profile = numpy.loadtxt(self.config['template_pulse_profile'],unpack=True)
		sys.stdout.write('[read_pulse_profile] %s is loaded.\n' % self.config['template_pulse_profile'])

		#shifts the main peak at phase 0.5
		main_peak_index  = self.template_profile.argmax()
		len_pulse        = len(self.template_profile)
		shift_pulse      = len_pulse-main_peak_index+len_pulse/2
		self.phase       = numpy.roll(self.phase,shift_pulse)
		self.template_profile = numpy.roll(self.template_profile,shift_pulse)
		if self.config['flag_template_error']:
			self.template_profile_error = numpy.roll(self.template_profile_error,shift_pulse)

		if debug:
			for i in range(number_of_sample_line):
				print "  -- ",
				print self.phase[i],
				print self.template_profile[i],
				if self.config['flag_template_error']:
					print self.template_profile_error[i]

	def run_fourier_transform(self):
		##takes the Furier transform
		sys.stdout.write('[run_fourier_transform] FFT transformed.\n')
		self.fourier_transformed = numpy.fft.fft(self.template_profile)
		if debug:
			for i in range(number_of_sample_line):
				print "  -- ",
				print self.fourier_transformed[i]
		self.length_of_FFT = len(self.fourier_transformed)
		if debug:
			sys.stdout.write('-- number of points of FFT: %d\n' % self.length_of_FFT)

	def remove_possoian_noise(self):
		##MAKE an FFT approximation of the pulse profile to remove Poissonian noise

		sys.stdout.write('[remove_possoian_noise] smoothing of the pulse profile.\n')		
		if self.config['flag_template_error']:
			if debug:
				sys.stdout.write('-- ** number_of_harmonics chisquare ** \n')
			for number_of_harmonics in range(1,self.length_of_FFT/2):
				##copies FFT into working array
				fft_copy = self.fourier_transformed.copy()
				##deletes harmonics higher than n_harm
				fft_copy[number_of_harmonics+1:self.length_of_FFT-1-number_of_harmonics]=0
 				##builds filterd profile
				filtered_profile = numpy.real(numpy.fft.ifft(fft_copy))
    			##computes sort of chi2
				chisquare = numpy.sum(((self.template_profile - filtered_profile)/self.template_profile_error)**2) / self.length_of_FFT	
				if debug: 
					sys.stdout.write('-- %d %.3e\n' % (number_of_harmonics, chisquare)) 
				##breaks loop when enough harmonics are included
				if chisquare <= self.config['chisquare_threshold']:
					break
			self.threshold_harmonics = number_of_harmonics
			sys.stdout.write('-- threshold_harmonics = %d\n' % self.threshold_harmonics)
			self.fft_copy = fft_copy
			self.filtered_profile = filtered_profile
		else:
			self.fft_copy = self.fourier_transformed.copy()
			self.filtered_profile = self.fft_copy

	def make_pulse_profile(self):
		sys.stdout.write('[make_pulse_profile] make pulse profle from FFT\n')
		##Make pulse profile with arbitrary number of bins (even) from FFT of the template
		
		#f1=np.zeros(n_bins,dtype='complex')
		#f1[0:len(ft)/2]=ft_f[0:len(ft)/2].copy()
		#f1[-len(ft)/2:]=ft_f[len(ft)/2:].copy()
		#f1*=float(n_bins)/len(r)
		#r1=np.real(np.fft.ifft(f1))
		#r1/=r1.mean()

		self.out_phase = numpy.linspace(0.,1.,self.config['number_of_bins'])
		self.out_fourier_transformed  = numpy.zeros(self.config['number_of_bins'],dtype='complex')
		self.out_fourier_transformed[0:self.length_of_FFT/2]=self.fft_copy[0:self.length_of_FFT/2].copy()
		self.out_fourier_transformed[-self.length_of_FFT/2:]=self.fft_copy[self.length_of_FFT/2:].copy()
		self.out_fourier_transformed *= float(self.config['number_of_bins'])/len(self.template_profile)
		self.out_pulse_profle = numpy.real(numpy.fft.ifft(self.out_fourier_transformed))
		#self.out_pulse_profle /= self.out_pulse_profle.mean()

	def simulate_pulse_random(self, exposure_ks):
		sys.stdout.write('[simulate_pulse_random] \n')
		# http://pacocat.com/?p=596

		exposure = exposure_ks * 1e+3
		self.config['number_of_pulsed_photons'] = int(exposure * self.config['source_rate'] * self.config['pulsed_fraction']) 
		self.config['number_of_DC_photons'] = int(exposure * (self.config['background_rate'] + self.config['source_rate'] * (1.-self.config['pulsed_fraction'])))
 	  	self.config['number_of_background_photons'] = int(exposure * self.config['background_rate'])
		self.config['number_of_photons'] = self.config['number_of_pulsed_photons'] + self.config['number_of_DC_photons']
		
		func = interp1d(self.out_phase,self.out_pulse_profle,kind='linear') # interpolate
		func_range = [min(self.out_phase),max(self.out_phase)] # holizontal range [x_min, x_max]
		func_max = max(self.out_pulse_profle) # vertical range		

		tstart_list = []
		tstop_list  = []
		gti_list    = []
		n = 0; exp = 0.; on_source = 0.0
		while exp <= exposure:
			tstart = n * self.config['orbit_period_min'] * 60.0 
			tstop  = tstart + self.config['orbit_period_min'] * 60.0 * self.config['gti_fraction']
			exp   += tstop - tstart 
			tstart_list.append(float(tstart))
			tstop_list.append(float(tstop))
			gti_list.append([tstart,tstop])
			n += 1 
		on_source = tstop 

		dump = '%d events (number of events) are generated...\n' % self.config['number_of_photons']
		sys.stdout.write(dump)

		number_of_cycle = int((on_source)/(self.config['period']*1e-3))
		sys.stdout.write('number of pulse cycle: %d (exp=%.3f ks, onSource=%.3f ks, P=%.3f ms)\n' % (
			number_of_cycle,
			exposure_ks, on_source/1e+3, 
			self.config['period']))

		sys.stdout.write("[%s]" % (" " * (TOOLBAR_WIDTH+1)))
		sys.stdout.flush()
		sys.stdout.write("\b" * (TOOLBAR_WIDTH+2)) # return to start of line, after '['
		TOOLBAR_BASE = int(self.config['number_of_photons'] / TOOLBAR_WIDTH)

		#phase_list = []
		#time_list  = []

		time_phase_list = []

		# pulsed signal generation 
		i = 0
		while i < self.config['number_of_pulsed_photons']: 
			phase  = randomgen(func,func_range,func_max)
			if self.config['time_noise_gauss_sigma'] > 0.0:
				phase += numpy.random.normal(0.0,self.config['time_noise_gauss_sigma']/self.config['period'])
			while phase < 0.0:
				phase += 1.0
			while phase > 1.0:
				phase -= 1.0 
			event_time = (float(numpy.random.randint(0,number_of_cycle))+phase) * self.config['period']*1e-3

			flag = False 
			for gti in gti_list:
				if event_time >= gti[0] and event_time <= gti[1]:
					flag = True
					break

			if flag == False:
				continue

			#phase_list.append(phase) 
			#time_list.append(event_time)
			time_phase_list.append([event_time,phase])
		 	i += 1

			if i % TOOLBAR_BASE == 0:
				sys.stdout.write("-")
		 		sys.stdout.flush()	
		sys.stdout.write('\n')

		# DC signal generation 
		i = 0
		while i < self.config['number_of_DC_photons']: 
			phase  = numpy.random.rand()
			gtinum = numpy.random.randint(len(gti_list))
			event_time = gti_list[gtinum][0] + phase * self.config['period']*1e-3

			#phase_list.append(phase) 
			#time_list.append(event_time)
			time_phase_list.append([event_time,phase])		
		 	i += 1

			if i % TOOLBAR_BASE == 0:
				sys.stdout.write("-")
		 		sys.stdout.flush()	
		sys.stdout.write('\n')

		time_phase_list.sort()
		time_list  = [ i[0] for i in time_phase_list]
		phase_list = [ i[1] for i in time_phase_list]

		if not os.path.exists(self.config['outdir_path']):
			os.system('mkdir -p %s' % self.config['outdir_path'])

		basename = '%s/%s_%dus_%dks' % (
			self.config['outdir_path'],
			self.config['outname'],
			1000.0*float(self.config['time_noise_gauss_sigma']),
			exposure_ks
			)
		qdpfile = '%s.qdp' % basename
		fitsfile = '%s.fits' % basename

		# ==========================
		# FITS 
		# ==========================
		hdu0 = pyfits.PrimaryHDU() 

		fits_col_TIME = pyfits.Column(name='TIME',format='1D',array=numpy.array(numpy.sort(time_list)),unit='s')
		fits_coldefs  = pyfits.ColDefs([fits_col_TIME])
		hdu1          = pyfits.new_table(fits_coldefs)
		hdu1.name     = 'EVENTS'

		fits_col_START = pyfits.Column(name='START',format='1D',array=numpy.array(tstart_list),unit='s')
		fits_col_STOP  = pyfits.Column(name='STOP',format='1D',array=numpy.array(tstop_list),unit='s')
		#fits_col_START = pyfits.Column(name='START',format='1D',array=numpy.array([0.0]),unit='s')
		#fits_col_STOP  = pyfits.Column(name='STOP',format='1D',array=numpy.array([on_source]),unit='s')		
		fits_coldefs2  = pyfits.ColDefs([fits_col_START,fits_col_STOP])
		hdu2           = pyfits.new_table(fits_coldefs2)
		hdu2.name     = 'GTI'		

		for header in [hdu0.header, hdu1.header, hdu2.header]:
			header.append(('TIMESYS','TT',''),end=True)
			header.append(('TIMEREF','LOCAL',''),end=True)
			header.append(('TIMEUNIT','s',''),end=True)
			header.append(('TSTART',0.0,''),end=True)
			header.append(('TSTOP',on_source,''),end=True)
			#header.append(('MJDREFI','51544','MJD reference day'),end=True)
			#header.append(('MJDREFF','7.428703703703700E-04','MJD reference (fraction of day)'),end=True)
			header.append(('TIMEDEL','1e-6','finest time resolution (time between frames'),end=True)
			header.append(('TIMEPIXR','0.0','0:times refer to the beginning of bin, 0.5:midd'),end=True)				
			header.append(('TIMEZERO','0.0','Time Zero'),end=True)				
		
		os.system('rm -f %s' % fitsfile)
		hdulist = pyfits.HDUList([hdu0, hdu1, hdu2])
		hdulist.writeto(fitsfile)
		#hdu = pyfits.PrimaryHDU(numpy.array(time_list))
		#
		#hdulist.writeto('test.fits')

		self.simulated_pulse = numpy.array(numpy.histogram(phase_list,bins=self.config['number_of_bins'],range=(0.0,1.0),normed=False, weights=None)[0])
		self.simulated_pulse_error = numpy.sqrt(self.simulated_pulse)
		fout = open(qdpfile,'w')
		fout.write('read serr 2\ncpd /xw\n')
		pcofile = '%s.pco' % (os.path.splitext(qdpfile)[0])
		fout.write('@%s\n' % pcofile)
		for i in range(self.config['number_of_bins']):
			dump = '%f %f %f \n' % (self.out_phase[i],
				self.simulated_pulse[i],
				self.simulated_pulse_error[i])
			fout.write(dump)
#		fout.write('no no no\n')
#		for i in range(self.config['number_of_bins']):
#			dump = '%f %f 0.0 \n' % (self.out_phase[i],
#				self.offseted_pulse_profile[i]*self.config['number_of_photons']/sum(self.offseted_pulse_profile))
#			fout.write(dump)		
		fout.close()		

		fout = open(pcofile,'w')
		title   = '%s ' % basename
		title  += '(P=%.3f s)' % (self.config['period']/1e+3)
		ftitle  = 'T=%.1f ks ' % exposure_ks
		ftitle += '(On=%.1f ks) ' % (on_source/1e+3)
		ftitle += 'pf=%d%% ' % (100.0*self.config['pulsed_fraction'])
		ftitle += 'R=%.2f ' % (self.config['source_rate'])
		ftitle += 'B=%.2f ' % (self.config['background_rate'])
		ftitle += 'jit.=%d \gms ' % (1000.0*float(self.config['time_noise_gauss_sigma']))
		ftitle += 'nbin=%d ' % (self.config['number_of_bins'])
#		ftitle += 'pulsed:%d cnts ' % (self.config['number_of_pulsed_photons'])
		dump  = 'skip on\nmark on\nmark 17 on 1\nerror off 2\nline on 2\nmark off 2\n'
		dump += 'la t %s\n' % title
		dump += 'la f\ntime off\nlab x Phase\nlab y Counts\n'
		dump += 'r x %.3f %.3f\n' % (-0.5/self.config['number_of_bins'],1+0.5/self.config['number_of_bins'])
		dump += 'r y 0 %f\n' % (1.2*self.simulated_pulse.max())
		#dump += 'r y 0 150\n' 
		#dump += 'mo gaus cons\n'
		#dump += '0.5\n'
		#dump += '0.01\n'
		#dump += '10\n'
		#dump += '20\n'
		#dump += 'fit\n'
		dump += 'line step on\n'
		dump += 'la f %s\n' % ftitle
		dump += 'la 1 P 0.0 %.1f LI 0 1 CO 4 LS 2 \n' % (self.config['number_of_background_photons']/self.config['number_of_bins'])
		dump += 'lwid 3 on 2'
		fout.write(dump)
		fout.close()

		psfile = '%s.ps' % basename
		cmd  = 'qdp %s <<EOF\n' % qdpfile
		cmd += 'hard %s/cps\n' % psfile
		cmd += 'quit\n'
		cmd += 'EOF\n'
		cmd += 'ps2pdf %s\n' % psfile
		cmd += 'mv %s.pdf %s\n' % (os.path.splitext(os.path.basename(psfile))[0],os.path.dirname(psfile))
		#cmd += 'open %s/%s.pdf \n' % (os.path.dirname(psfile),os.path.splitext(os.path.basename(psfile))[0])
		print cmd; os.system(cmd)

		f = open('%s.yaml' % (basename),'w')
		f.write(yaml.dump(self.config))
		f.close()

	def run(self, number_of_photons):
		self.show_parameters()
		self.read_pulse_profile()
		self.run_fourier_transform()
		self.remove_possoian_noise()
		self.make_pulse_profile()
		self.simulate_pulse_random(number_of_photons)

if __name__ == "__main__":

	if len(sys.argv) != 3:
		sys.stderr.write('%s input.yaml exposure_ks\n' % sys.argv[0])
		quit()
	yamlfile    = sys.argv[1]
	exposure_ks = sys.argv[2]
	pulse = PulseProfile(yamlfile)
	pulse.run(float(exposure_ks))

