####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for producing approximate ME results for e+e- -> qqBar with massless particles."""

##Import required modules:##
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '', 'Tools'))
import math
import assertions
import precision
import constants
import counters
import dataLoggers
import particleData
import fourVectors
import particles
import kinematics
import quarkPairs
import LHEFHandlers

print "\n///////////////////"
print "Loading MEs module:"
print "///////////////////\n"

##Set up counters:##
quarkCounters = []
for i in particleData.knownParticles.get_known_quarks():
	quarkCounters.append(counters.counter(1))
createdThetas = dataLoggers.dataLogger()
createdQuarkCodes = dataLoggers.dataLogger()

##Functions:##

def check_is_qqBarMEGenerator(toCheck):
	"""A function to check for an instance of the qqBarMEGenerator class."""
	return isinstance(toCheck,qqBarMEGenerator)

def produce_electron_pair(S123):
	"""A function to produce an electron and positron in opposite directions."""
	__cE = constants.cut_off_energy()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	##Define beams to be along z: works well for Pythia 8.2.
	__E = math.sqrt(S123)/2.0
	__pV1 = fourVectors.fourVector(__E,0.0,0.0,__E)
	__pV2 = fourVectors.fourVector(__E,0.0,0.0,-__E)
	__eCode = particleData.knownParticles.get_code_from_name("electron")
	__p1 = particles.particle(__eCode,__pV1,[0,0],[3,4],[0,0],-1)
	__p2 = particles.particle(-__eCode,__pV2,[0,0],[3,4],[0,0],-1)
	__p1.set_unique_ID(1)
	__p2.set_unique_ID(2)
	return [__p1,__p2]

def produce_quark_directions(S123,theta):
	"""A function to produce a set of opposite directional quark vectors in opposite directions."""
	__cE = constants.cut_off_energy()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	assert((type(theta) == float) and (0.0 <= theta) and (theta <= math.pi))
	__phi = kinematics.get_random_phi()
	__x = math.sin(theta)*math.cos(__phi)
	__y = math.sin(theta)*math.sin(__phi)
	__z = math.cos(theta)
	__direction = fourVectors.fourVector(0.0,__x,__y,__z)
	__oppositeDirection = __direction.copy()
	__oppositeDirection *= -1.0
	return __direction, __oppositeDirection

def scaled_quark_directions(S123,theta):
	"""A function to return two quark directional vectors scaled to match a given S123 value."""
	__cE = constants.cut_off_energy()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	assert((type(theta) == float) and (0.0 <= theta) and (theta <= math.pi))
	__E = math.sqrt(S123)/2.0
	__dV1, __dV2 = produce_quark_directions(S123,theta)
	__magnitude = __dV1.calculate_cartesian_magnitude()
	__scale = __E/__magnitude
	__dV1 *= __scale
	__dV2 *= __scale
	__dV1[0], __dV2[0] = __E, __E
	##Check massless:
	assert (precision.check_numbers_equal(__dV1*__dV1,0.0) and precision.check_numbers_equal(__dV2*__dV2,0.0))
	return __dV1, __dV2

def produce_quark_pair(S123,posQuarkCodes):
	"""A function to produce a particle and anti-particle in opposite directions."""
	__cE = constants.cut_off_energy()
	__knownQuarks = particleData.knownParticles.get_known_quarks()
	assert ((type(S123) == float) and (S123 > __cE*__cE))
	assert (type(posQuarkCodes) == list)
	for code in posQuarkCodes:
		assert (type(code) == int)
		assert (code in __knownQuarks)
	__qCode, __qTheta = quarkPairs.get_code_and_theta(S123,posQuarkCodes)
	quarkCounters[__qCode - 1].count() ##-1 for list index.
	createdThetas.store(__qTheta)
	createdQuarkCodes.store(__qCode)
	__pV1, __pV2 = scaled_quark_directions(S123,__qTheta)
	__p1 = particles.particle(__qCode,__pV1,[1,2],[0,0],[501,0],1)
	__p2 = particles.particle(-__qCode,__pV2,[1,2],[0,0],[0,501],1)
	__p1.set_unique_ID(3)
	__p2.set_unique_ID(4)
	__p1.set_produced_at(S123)
	__p2.set_produced_at(S123)
	return [__p1,__p2]

##Classes:##

class qqBarMEGenerator(object):
	"""A class for generating approximate ME results for e+e- -> qqBar with massless particles."""

	def __init__(self,S123,quarkCodes,numEvents):
		"""A function to initalise an approximate ME generator for massless particles."""
		##Only expects quark codes not anti-quark codes at these are implicit.
		assert ((type(S123) == float) and (S123 > 0))
		assert (type(quarkCodes) == list)
		__knownQuarks = particleData.knownParticles.get_known_quarks()
		self.__electronCode = particleData.knownParticles.get_code_from_name('electron')
		for code in quarkCodes:
			assert (type(code) == int)
			assert (abs(code) in __knownQuarks)
		assert ((type(numEvents) == int) and (numEvents > 0))
		##list(set()) removes any duplicates in case there were any.
		self.__S123, self.__quarkCodes, self.__numEvents = S123, list(set(quarkCodes)), numEvents

	def run(self):
		"""A function to run the qqBar ME generator approximation."""
		##Initialise writer
		__writer = LHEFHandlers.LHEFMEWriter('qqBar',self.__numEvents,self.__S123,self.__quarkCodes)
		__trials, __muf2, __mur2 = 1, self.__S123, self.__S123
		__processID, __weight = 0, 1
		__scale, __alphaEM = 0.0,0.0
		__alphaS = 0.0
		##Generate and output events.
		__i, __printEvery = 0, self.__numEvents/10
		if __printEvery > 1000: ##Slow enough to want to see something is happening!
			__printEvery = 1000
		elif (__printEvery == 0): ##i.e < 10 events.
			__printEvery = 1
		while (__i < self.__numEvents):
			__initial = produce_electron_pair(self.__S123)
			__final = produce_quark_pair(self.__S123,self.__quarkCodes)
			__writer.add_event(__trials,__muf2,__mur2,__processID,__weight,__scale,__alphaEM,__alphaS,__initial,__final)
			if (__i%__printEvery == 0):
				if (__i == 0):
					print "\nStart generating..."
				else:
					print "Reached event", __i
			__i += 1 ##Iterate
		##End output.
		print "\nSaving..."
		__writer.save()
		assertions.pause(__name__)
		print "\n---------------------"
		print "Quark content report:"
		print "---------------------"
		__sumPercent = 0.0
		for __qIndex, __counter in enumerate(quarkCounters):
			__qName = particleData.knownParticles.get_name_from_code(__qIndex + 1)
			__qPercent = __counter.counted()*100.0/self.__numEvents
			__sumPercent += __qPercent
			print str(__qPercent) + "% are " + __qName + "s."
		print "This adds up to " + str(__sumPercent) + "%."

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##

	testMode = None
	while (type(testMode) != str):
		testMode = str(raw_input("\nEnter y to start test mode. Any other input will start the generator: "))
	
	##Setup shared data here:##
	__Zcode = particleData.knownParticles.get_code_from_name('Z-boson')
	__Mz = particleData.knownParticles.get_mass_from_code(__Zcode)
	
	if not ((testMode == 'y') or (testMode == 'Y')):
		##Import additional required modules here:##

		##Begin running:##
		print "\n----------------------------------------------------------------------"
		print "----------------------------------------------------------------------\n"
		print "///////////////////"
		print "Running MEs module:"
		print "///////////////////"
		assertions.pause(__name__)

		##Setup generator here:##
		__S123, __qCodes = None, None
		while ((type(__S123) != float) and (type(__qCodes) != list)):
			try:
				__numEvents = int(raw_input("\nEnter the number of events to generate: "))
				__S123 = raw_input("\nEnter S123 to use or leave blank for default: ")
				if (__S123 == ""):
					__S123 = __Mz*__Mz
				else:
					__S123 = float(__S123)
				__qCodes = raw_input("\nEnter quark codes to use seperated by a space or leave blank for default: ")
				if (__qCodes == ""):
					__qCodes = constants.active_q_codes()
				else:
					__qCodes = [int(__x) for __x in __qCodes if __x != ' ']
			except:
				pass
		##Run##
		theqqBarMEGenerator = qqBarMEGenerator(__S123,__qCodes,__numEvents)
		theqqBarMEGenerator.run()

		###Plot produced thetas:##
		print "\n--------------------------------------------------\n"
		__yN = raw_input("Enter 'y' to plot produced theta values or leave blank to continue: ")
		if ((__yN == "y") or (__yN == "Y")):
			##Import required modules for plotting:##
			import matplotlib.ticker as mtick
			from matplotlib import pyplot
			import numpy

			##Setup required test values:##
			__numberThetas = 100000
			__cosThetaRange = [math.cos(i*math.pi/(__numberThetas - 1.0)) for i in range(0,__numberThetas)]
			__numCosThetaIts = float(__numberThetas)
			__numberPerBin = 100
			__numberBins = int(__numCosThetaIts)/__numberPerBin
			__testQNames = [r"$\rm{d-quark}$",r"$\rm{u-quark}$",r"$\rm{s-quark}$",r"$\rm{c-quark}$",r"$\rm{b-quark}$",r"$\rm{t-quark}$"]

			##Prepare and plot graphs:##
			print "\nTesting produced theta values:\n"
			__cosThetaProbs = []
			__allCosThetaHists, __allCosThetaBins, __allBinCentres = [], [], []
			__cosThetasOut = []
			for __aQCode in __qCodes:
				print "Using code " + str(__aQCode) +":"
				__currentCosThetaProbs = [quarkPairs.calc_f_for_cos_theta(__S123,__aQCode,__aCosTheta) for __aCosTheta in __cosThetaRange]
				thetaProbsSum = sum(__currentCosThetaProbs)
				__cosThetaProbs.append([__i/float(thetaProbsSum) for __i in __currentCosThetaProbs])
				print "Normalised expected values sum to:", sum(__cosThetaProbs[-1])
				__loggedCodes = createdQuarkCodes.output()
				__loggedThetas = createdThetas.output()
				##Split out the theta values corresponding to this quark type:
				__currentCosThetas = [math.cos(__loggedThetas[__i]) for __i, __v in enumerate(__loggedCodes) if __v == __aQCode]
				__cosThetasOut.append(__currentCosThetas)
				__bins = numpy.linspace(-1.0,1.0,__numberBins+1)
				__cosThetasHists, __cosThetaBins = numpy.histogram(__cosThetasOut[-1],__bins)
				##Normalise to 1. Have to divide by the number per bin for plotting the average.
				__allCosThetaHists.append([__aNumberInBin/float(len(__currentCosThetas)*__numberPerBin) for __aNumberInBin in __cosThetasHists])
				print "Normalised generated values sum to:", sum(__allCosThetaHists[-1])*__numberPerBin ##As divided each by above.
				__allCosThetaBins.append(__cosThetaBins)
				__allBinCentres.append([])
				for __i in range(len(__allCosThetaBins[-1])-1): ##-1 for number of centres
					__allBinCentres[-1].append((__allCosThetaBins[-1][__i] + __allCosThetaBins[-1][__i+1])/2.0)
			__xaxes = [r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$",r"$\rm{cos\ \theta}$"]
			__yaxes = [r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$"]
			__yaxes += [r"$\rm{P}(\rm{cos\ \theta})$",r"$\rm{P}(\rm{cos\ \theta})$"]
			__yMaxs = [0.000025,0.000025,0.000025,0.000025,0.000025,0.000025]
			__supTitle = r"$\rm{Monte\ Carlo\ sampling\ of\ cos\ \theta\ for}\ e^{+}e^{-}\ \rightarrow\ q\bar{q}\ \rm{using\ "
			__supTitle += str(int(__numCosThetaIts)) + r"\ iterations}$"
			__figure,__axes = pyplot.subplots(3,2)
			__axes = __axes.ravel()
			for __idx,__ax in enumerate(__axes):
				if not ((__idx + 1) in __qCodes):
					continue
				__ax.plot(__allBinCentres[__idx], __allCosThetaHists[__idx], linestyle = "solid", color = "blue", linewidth = 2)
				__ax.plot(__cosThetaRange, __cosThetaProbs[__idx], linestyle = "solid", color = "red", linewidth = 2)
				__ax.set_title(__testQNames[__idx])
				__ax.set_xlabel(__xaxes[__idx])
				__ax.set_ylabel(__yaxes[__idx])
				__ax.axis([-1,1,0,__yMaxs[__idx]])
				__ax.set_yticks([0.0,0.00001,0.00002,0.00003],minor=False)
				__ax.yaxis.set_major_formatter(mtick.FixedFormatter([r"$0.0$",r"$1.0$",r"$2.0$",r"$\times\ 10^{-3}$"]))
				__ax.set_xticks([-1.0,-0.5,0.0,0.5,1.0],minor=False)
				__ax.xaxis.set_major_formatter(mtick.FixedFormatter([r"$\minus1.0$",r"$\minus0.5$",r"$0.0$",r"$0.5$",r"$1.0$"]))
			pyplot.suptitle(__supTitle,fontsize = "16")
			__line1 = pyplot.Line2D((0,1),(0,0), color="blue", linewidth = 2)
			__line2 = pyplot.Line2D((0,1),(0,0), color="red", linewidth = 2)
			__lines, __figStrs = [__line1,__line2], [r"$\rm{Generated}$",r"$\rm{Expected}$"]
			__figure.legend(__lines, __figStrs, bbox_to_anchor=[0.5, 0.05],loc='center', ncol=2)
			pyplot.tight_layout()
			##Space main title out to prevent overlapping and allow space for legend below:
			pyplot.subplots_adjust(top=0.85,bottom=0.17)
			assertions.show_graph()
			print "\nFinished testing produced theta values."
			assertions.pause(__name__)

		##Done running:##
		print "\n---------------------------------------------\n"
		print "////////////////////////////"
		print "Finished running MEs module!"
		print "////////////////////////////"

	else:
		##Import additional required modules here:##
		from matplotlib import pyplot
		import numpy

		##Begin testing:##
		print "\n----------------------------------------------------------------------"
		print "----------------------------------------------------------------------\n"
		print "///////////////////"
		print "Testing MEs module:"
		print "///////////////////"
		assertions.pause(__name__)

		##Setup test code here:##

		##Done testing:##
		print "\n---------------------------------------------\n"
		print "/////////////////////////////"
		print "Finished checking MEs module!"
		print "/////////////////////////////"