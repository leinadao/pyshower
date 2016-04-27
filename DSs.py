####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling the dipole showering of the events in a LHEF XML file."""

##Import required modules:##
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '', 'Tools'))
import assertions
import precision
import counters
import particleData
import kinematics
import sudakovs
import chains
import showers
import LHEFHandlers

print "\n///////////////////"
print "Loading DSs module:"
print "///////////////////\n"

##Functions:##

def check_is_controller(toCheck):
	"""A function to check for an instance of the controller class."""
	return isinstance(toCheck,controller)

##Classes:##

class controller(object):
	"""A class for handling a dipole shower."""

	def __init__(self,fileName,activeQCodes):
		"""A function to initiate a controller for showering an LHEF XML file."""
		assert (type(fileName) == str)
		assert (type(activeQCodes) == list)
		for qCode in activeQCodes:
			assert qCode in particleData.knownParticles.get_known_quarks()
		self.__activeQCodes = activeQCodes
		if not (fileName[-4:] == '.lhe'):
			self.__toLoad = fileName + '.lhe'
		else:
			self.__toLoad = fileName

	def run(self):
		"""A function to run a shower for the given LHEF file."""
		##For e+e- -> qqBar (massless).
		__reader = LHEFHandlers.LHEFReader(self.__toLoad)
		self.__numEvents = __reader.get_number_events()
		__writer = LHEFHandlers.LHEFShowerWriter(__reader,self.__numEvents)
		__alertEvery = self.__numEvents/10
		if __alertEvery > 1000: ##Slow enough to want to see something is happening!
			__alertEvery = 1000
		elif (__alertEvery == 0): ##i.e < 10 events.
			__alertEvery = 1
		for __eventIndex in range(self.__numEvents):
			if (__eventIndex%__alertEvery == 0):
				if (__eventIndex == 0):
					print "Begin showering...\n"
				else:
					print "Showering event", __eventIndex
			__MEParticle1 = __reader.get_event_particle(__eventIndex,0)
			__MEParticle2 = __reader.get_event_particle(__eventIndex,1)
			__EIn = __MEParticle1[0] + __MEParticle2[0]
			__S123In = kinematics.Sijk([__MEParticle1.get_four_momentum(), __MEParticle2.get_four_momentum()])
			__particlesIn = [__MEParticle1,__MEParticle2]
			__showerParticle1 = __reader.get_event_particle(__eventIndex,2)
			__showerParticle2 = __reader.get_event_particle(__eventIndex,3)
			__showeri = showers.qqBarShower(__showerParticle1,__showerParticle2,self.__activeQCodes)
			__showeri.run_shower()
			__particlesOut = __showeri.export_results()
			__EOut = 0.0
			__fourVectorsOut = []
			for __p in __particlesOut:
				__ePV = __p.get_four_momentum()
				__EOut += __ePV[0]
				__fourVectorsOut.append(__ePV)
			assert precision.check_numbers_equal(kinematics.Sijk(__fourVectorsOut),__S123In)
			assert precision.check_numbers_equal(__EIn,__EOut)
			__trials, __muf2, __mur2 = __reader.get_event_trials(__eventIndex), __reader.get_event_muf2(__eventIndex), __reader.get_event_mur2(__eventIndex)
			__processID, __weight = __reader.get_event_process_ID(__eventIndex), __reader.get_event_weight(__eventIndex)
			__scale, __alphaEM = __reader.get_event_scale(__eventIndex), __reader.get_event_alphaEM(__eventIndex)
			__alphaS = __reader.get_event_alphaS(__eventIndex)
			__writer.add_event(__trials,__muf2,__mur2,__processID,__weight,__scale,__alphaEM,__alphaS,__particlesIn,__particlesOut)
		__writer.save()
		__numGluonsSplit = counters.gluonSplitCounter.counted()
		print "\nThere were", counters.gluonProdCounter.counted(), "gluons produced!"
		print "\nThere were", __numGluonsSplit, "gluons split!"
		print "\nThere were", counters.photonProdCounter.counted(), "photons produced!"
		print "\nThere were", counters.kPerpProdWarningCounter.counted(), "warnings for E1 or E3 < E2!\n"
		print "\n---------------------"
		print "Quark content report:"
		print "---------------------"
		__sumPercent = 0.0
		for __aQCode in [1,2,3,4,5,6]:
			__qName = particleData.knownParticles.get_name_from_code(__aQCode)
			if precision.check_numbers_equal(__numGluonsSplit,0.0):
				__qPercent = 0.0
			else:
				__qPercent = chains.producedQuarkCodes.output().count(__aQCode)*100.0/__numGluonsSplit
			__sumPercent += __qPercent
			print str(__qPercent) + "% are " + __qName + "s."
		print "This adds up to " + str(__sumPercent) + "%."
		assertions.pause(__name__)

##Module running code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import pyperclip
	from matplotlib import pyplot
	import constants

	##Begin running:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////"
	print "Running DSs module:"
	print "///////////////////"
	assertions.pause(__name__)
	
	def plot_energy_ratios():
		"""A function to plot the energy ratios generated in the kinematics module throughout."""
		E1s = kinematics.e1s.output()
		E2s = kinematics.e2s.output()
		E3s = kinematics.e3s.output()
		ratio1s, ratio2s = [], []
		for i, x in enumerate(E1s):
			ratio1s.append([E2s[i]/E1s[i]])
			ratio2s.append([E2s[i]/E3s[i]])
		ratios1s = [x[0] for (y,x) in sorted(zip(E2s,ratio1s), key=lambda pair: pair[0])]
		ratios2s = [x[0] for (y,x) in sorted(zip(E2s,ratio2s), key=lambda pair: pair[0])]
		E2s.sort()
		pyplot.figure()
		pyplot.title(r"$The\ energy\ ratios\ produced\ throughout$")
		pyplot.xlabel(r"$E_{2}$")
		pyplot.ylabel(r"$Ratio E_{2}/E_{i}$")
		pyplot.yscale('log')
		pyplot.scatter(E2s,ratio1s,linewidth = 2, label = r"$Ratio E_{2}/E_{1}$")
		pyplot.scatter(E2s,ratio2s,linewidth = 2, label = r"$Ratio E_{2}/E_{3}$")
		pyplot.legend()
		assertions.show_graph()

	def plot_xs():
		"""A function to plot x1 vs x3 generated in the kinematics module throughout."""
		x1s = kinematics.tX1s.output()
		x3s = kinematics.tX3s.output()
		pyplot.figure()
		pyplot.xlabel(r"$x_{1}$",fontsize = 22)
		pyplot.ylabel(r"$x_{3}$",fontsize = 22)
		pyplot.xlim(-0.05,1.05)
		pyplot.ylim(-0.05,1.05)
		pyplot.xticks(fontsize = 15)
		pyplot.yticks(fontsize = 15)
		pyplot.scatter(x1s,x3s)
		assertions.show_graph()

	##Setup here:##
	__thePath, __activeQCodes = None, None
	while ((type(__thePath) != str) or (__thePath == "")):
		try:
			__thePath = str(raw_input("\nEnter the relative path of the LHEF input or leave blank to use clipboard: "))
			if __thePath == "":
				__thePath = str(pyperclip.paste())
			__activeQCodes = raw_input("\nEnter quark codes to use seperated by a space or leave blank for default: ")
			if (__activeQCodes == ""):
				__activeQCodes = constants.active_q_codes()
			else:
				__activeQCodes = [int(__x) for __x in __activeQCodes if __x != ' ']
		except:
			pass
	##Run##
	theController = controller(__thePath,__activeQCodes)
	theController.run()
	__crossSecsGluSplit = sudakovs.crossSecsGluSplit.output()
	__crossSecsGluProd = sudakovs.crossSecsGluProd.output()
	__ratios = [__crossSecsGluSplit[i]/__crossSecsGluProd[i] for i,x in enumerate(__crossSecsGluProd)]
	__averageRatio = sum(__ratios)/float(len(__ratios))
	print "\nThe average ratio of the gluon splitting cross section to the gluon production cross section was:", __averageRatio, "\n"
	plot_energy_ratios()
	plot_xs()

	##Done:##
	print "\n---------------------------------------------\n"
	print "////////////////////////////"
	print "Finished running DSs module!"
	print "////////////////////////////"