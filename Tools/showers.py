####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling dipole showering."""

##Import required modules:##
import assertions
import counters
import particleData
import particles
import dipoles
import chains
import resultsContainers

print "\n///////////////////////"
print "Loading showers module:"
print "///////////////////////\n"

##Functions:##

def check_is_qqBarShower(toCheck):
	"""A function to check for an instance of the qqBarShower class."""
	return isinstance(toCheck,qqBarShower)

##Classes:##

class qqBarShower(object):
	"""A class for handling a dipole shower."""

	def __init__(self,particle1,particle2,activeQCodes,history=False):
		"""A function to initiate a qqBar shower."""
		assert particles.check_is_particle(particle1)
		assert particles.check_is_particle(particle2)
		assert (particle1.get_code() == -1*particle2.get_code())
		assert (type(activeQCodes) == list)
		for qCode in activeQCodes:
			assert qCode in particleData.knownParticles.get_known_quarks()
		self.__activeQCodes = activeQCodes
		##Provide unique particle identifiers for produced particles per event.
		self.__particleCounter = counters.counter(5) ##Start at 4 as e+,e-,q,qBar will be 1,2,3,4 for qqBar shower.
		self.__colourCounter = counters.counter(502) ##Start at 501 to differentiate but 501 taken already by qqBar.
		##Prepare particles etc.
		self.__startParticle1 = particle1.copy()
		self.__startParticle2 = particle2.copy()
		self.__startDipole = dipoles.dipole(self.__startParticle1.copy(),self.__startParticle2.copy())
		self.__startChain = chains.chain([self.__startParticle1.copy(),self.__startParticle2.copy()])
		self.__container = resultsContainers.resultsContainer()
		self.__isRun = False
		self.__showerList = [chains.chain([self.__startParticle1,self.__startParticle2])]
		self.__photonList = []
		self.__showerHistory = []

	def get_start_dipole():
		"""A function to return the dipole a shower began from."""
		return self.__startDipole

	def get_history():
		"""A function to return the history of a particle."""
		return self.__showerHistory

	def is_run(self):
		"""A function to determine whether the shower has been run."""
		return self.__isRun

	def codes_produced(self):
		"""A function to return the codes of all resulting particles."""
		assert self.__isRun
		__codesList = []
		for __chain in self.__showerList:
			for __dipole in __chain.get_chain_list():
				__codesList.append(__dipole[0].get_code())
			if not __chain.check_is_closed(): ##Get the extra particle on the RHS of the end dipole.
				__codesList.append(__chain[-1][1].get_code())
	 	return __codesList

	def number_particles_out(self,code = "all"):
		"""A function to count the number of resulting particles of a given type."""
		assert (code == "all" or (type(code) == int))
		assert self.__isRun
		__codesList = self.codes_produced()
		if (code == "all"):
			__number = len(__codesList)
		else:
			__number = __codesList.count(code)
	 	return __number

	def save(self):
		"""A function to save all of the resulting particles in a results container."""
		assert self.__isRun
		for __chain in self.__showerList:
			__chainList = __chain.get_chain_list()
			for __dipole in __chainList:
				self.__container.store(__dipole[0].copy())
			if not __chain.check_is_closed(): ##Already included at [0][0] if closed chain.
				self.__container.store(__chainList[-1][1].copy())
		for __photon in self.__photonList:
			self.__container.store(__photon.copy())

	def __getitem__(self, index):
		"""A function to get a chain by index after the shower has run."""
		assert (type(index) == int)
		assert self.__isRun
		return self.__showerList[index]

	def export_results(self):
		"""A function to export the list of resulting particles from the results container."""
		assert self.__isRun
		return self.__container.get_all()

	def all_chains_showered(self):
		"""A function to check if all chains are showered to completion."""
		for __chain in self.__showerList:
			if (not __chain.showering_completed()):
				return False
		return True

	def run_shower(self):
		"""A function to run the shower to completion."""
		#Currently shower history can't be turned on.
		while (not self.all_chains_showered()):
			__breakInnerLoops = False
			__numChains = len(self.__showerList)
			for __showerListIndex in range(__numChains): ##Iterating length not list to prevent continuing when list updated.
				__chain = self.__showerList[__showerListIndex]
				if (not __chain.showering_completed() and (not __breakInnerLoops)):
					__results = [4,None] ##Used to start while loop.
					while ((__results[0] != 0) and (not __breakInnerLoops)):
						__results = self.__showerList[__showerListIndex].evolve(self.__activeQCodes,self.__particleCounter,self.__colourCounter)
						if (__results[0] == 1): ##Gluon emission occured and no action required here.
							counters.gluonProdCounter.count()
						if (__results[0] == 2): ##Gluon splitting occured.
							counters.gluonSplitCounter.count()
							__splitBeforeIndex, __newMaxPperpSquared = __results[1]
							__chainPart1 = __chain.get_chain_list()[:(__splitBeforeIndex)] ##The 'middle' dipole is not needed.
							__chainPart2 = __chain.get_chain_list()[(__splitBeforeIndex + 1):]
							__orderedList1, __orderedList2 = [], []
							for __dipole in __chainPart1:
								__orderedList1.append(__dipole[0])
							__orderedList1.append(__chainPart1[-1][1]) ##Don't forget last particle.
							for __dipole in __chainPart2:
								__orderedList2.append(__dipole[0])
							__orderedList2.append(__chainPart2[-1][1]) ##Don't forget last particle.
							#Splitting of loops not yet treated here:
							__newChain1 = chains.chain(__orderedList1,False,__newMaxPperpSquared)
							__newChain2 = chains.chain(__orderedList2,False,__newMaxPperpSquared)
							self.__showerList[__showerListIndex] = __newChain1
							self.__showerList.insert(__showerListIndex + 1,__newChain2)
							__breakInnerLoops = True
						if (__results[0] == 3): ##Photon emission occured.
							counters.photonProdCounter.count()
							__producedPhoton = __results[1]
							self.__photonList.append(__producedPhoton.copy())
		self.__isRun = True
		self.save()

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import LHEFHandlers

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////"
	print "Testing showers module:"
	print "///////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	print "For qqBar dipole initiated shower:"
	testPath, eventIndex, particleIndex = None, None, None
	while ((type(testPath) != str) or (type(eventIndex) != int)):
		try:
			testPath = str(raw_input("\nEnter the relative path of the LHEF input or leave blank for default: "))
			if testPath == "":
				testPath = "../Events/belle.lhe"
			eventIndex = int(raw_input("\nEnter the index of the event in the LHEF file to use: "))
		except:
			pass
	testReader = LHEFHandlers.LHEFReader(testPath)
	testPDGID1 = testReader.get_event_particle_ID(eventIndex,2) ##First q, as 1 and 2 are e^+ and e^-.
	testPDGID2 = testReader.get_event_particle_ID(eventIndex,3)
	testVector1 = testReader.get_event_particle_four_momentum(eventIndex,2)
	testVector2 = testReader.get_event_particle_four_momentum(eventIndex,3)
	testParticle1 = particles.particle(testPDGID1,testVector1)
	testParticle2 = particles.particle(testPDGID2,testVector2)
	testShower1 = qqBarShower(testParticle1,testParticle2)

	##Test check_is_qqBarShower:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_qqBarShower:\n"
	print "Calling check_is_qqBarShower on instance: " , check_is_qqBarShower(testShower1)
	print "Calling check_is_qqBarShower on wrong instance: " , check_is_qqBarShower(testParticle1)
	print "Calling check_is_qqBarShower on 1.055: " , check_is_qqBarShower(1.055)
	print "Calling check_is_qqBarShower on 'word': " , check_is_qqBarShower('word')
	testResults = [check_is_qqBarShower(testShower1),check_is_qqBarShower(testParticle1)]
	testResults += [check_is_qqBarShower(1.055),check_is_qqBarShower('word')]
	if (sum(testResults) == 1):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_qqBarShower."
	assertions.pause(__name__)

	##Test shower class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////"
	print "Testing shower class:"
	print "/////////////////////"
	assertions.pause(__name__)

	##Test run_shower():##
	print "\n--------------------------------------------------\n"
	print "Testing run_shower():\n"
	testShower1.run_shower()
	print "\nFinished testing run_shower()."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////"
	print "Finished checking showers module!"
	print "/////////////////////////////////"