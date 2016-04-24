####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for storing results produced during dipole showering."""

##Import required modules:##
import particles

print "\n/////////////////////////////////"
print "Loading resultsContainers module:"
print "/////////////////////////////////\n"

##Functions:##

def check_is_resultsContainer(toCheck):
	"""A function to check for an instance of the resultsContainer class."""
	return isinstance(toCheck,resultsContainer)

##Classes:##

class resultsContainer(object):
	"""A class for storing results from dipole showering."""

	def __init__(self):
		self.__containerList = []

	def __getitem__(self,index):
		"""A function to get a result by index from a container."""
		return self.__containerList[index]

	def get_all(self):
		"""A function to get the list of all results in a container."""
		return self.__containerList

	def store(self,particleIn):
		"""A function to store a particle in a results container."""
		assert particles.check_is_particle(particleIn)
		self.__containerList.append(particleIn)

	def copy(self):
		"""A function to copy a results container."""
		__aCopy = resultsContainer()
		for __result in self.__containerList:
			__aCopy.store(__result.copy())
		return __aCopy

	def __eq__(self,otherRC):
		"""A function to check if two results containers are identical."""
		assert check_is_resultsContainer(otherRC)
		if (len(self.__containerList) != len(otherRC.get_all())): ##Quicker test if they aren't.
			return False
		elif ((len(self.__containerList) == 0) and (len(otherRC.get_all()) == 0)):
			return True ##Can only be equal if completely empty!
		for i, particle in enumerate(self.__containerList):
			if particle != otherRC[i]:
				return False
		return True

	def __ne__(self,otherRC):
		"""A function to check if two results containers are not identical."""
		assert check_is_resultsContainer(otherRC)
		return (not self.__eq__(otherRC))

	def __iadd__(self,toAdd): ##Used for += on results containers.
		"""A function for adding another results container to itself."""
		assert check_is_resultsContainer(toAdd)
		##The chances of all the four-vectors in genuinely different RC's being exactly equal is so low:
		assert self.__ne__(toAdd) ##Stops accidental duplication of results.
		__addList = toAdd.get_all()
		__combined = self.copy()
		for __i in range(len(__addList)):
			__combined.store(__addList[__i])
		return self

	def __add__(rCA, rCB): ##Used for a + b on results containers.
		"""A function for addition of results container objects."""
		assert (check_is_resultsContainer(rCA) and check_is_resultsContainer(rCB))
		##The chances of all the four-vectors in genuinetly different RC's being exactly equal is so low:
		assert rCA.__ne__(rCB) ##Stops accidental duplication of results.
		__combined = rCA.copy()
		__combined += rCB
		return __combined

	def __str__(self):
		"""A function to return a string of a resultsContainer object."""
		if (len(self.__containerList) == 0):
			return "[RC-[!EMPTY!]-]"
		__stringOfRC = "\n[RC-"
		for particle in self.__containerList:
			__stringOfRC += (str(particle) + " |\n")
		__stringOfRC += ("\n-]")
		return __stringOfRC

	def simple_str(self):
		"""A function to return a simplified string of a resultsContainer object."""
		if (len(self.__containerList) == 0):
			return "[RC-[!EMPTY!]-]"
		__simpleStringOfRC = "[RC-"
		for particle in self.__containerList:
			__simpleStringOfRC += particle.simple_str()
		__simpleStringOfRC += "]"
		return __simpleStringOfRC

	def __repr__(self):
		"""A function to return a string of a resultsContainer object."""
		if (len(self.__containerList) == 0):
			return "<resultsContainer([!EMPTY!])>"
		__reprOfRC = "\n<resultsContainer(["
		for particle in self.__containerList:
			__reprOfRC += (repr(particle) + " |\n")
		__reprOfRC += ("])>")
		return __reprOfRC

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import random
	import math
	import assertions
	import counters
	import fourVectors
	import particles

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////////////////"
	print "Testing resultsContainers module:"
	print "/////////////////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	##Generate random test four-vectors:##
	testFourVectors = {'a':None,'b':None,'c':None,'d':None,'e':None}
	##Prevent the random vectors being the same as can't boost into frame of massless particle.
	testVector1, testVector2 = fourVectors.fourVector(0,0,0,0), fourVectors.fourVector(0,0,0,0)
	selectionFinished = False
	while (not selectionFinished):
		for testFourVector in testFourVectors:
			x1 = random.randrange(1,100)
			x2 = random.randrange(0,100)
			x3 = random.randrange(0,100)
			##Make them massless.
			x0 = math.sqrt((x1*x1) + (x2*x2) +(x3*x3))
			testFourVectors[testFourVector] = fourVectors.fourVector(x0,x1,x2,x3)
		testVector1 = testFourVectors['a'].copy()
		testVector2 = testFourVectors['b'].copy()
		testVector3 = testFourVectors['c'].copy()
		testVector4 = testFourVectors['d'].copy()
		testVector5 = testFourVectors['e'].copy()
		##Want five independent energy-momentum four-vectors.
		if (not fourVectors.check_different_direction(testVector1,testVector2)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector1,testVector3)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector1,testVector4)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector1,testVector5)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector2,testVector3)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector2,testVector4)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector2,testVector5)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector3,testVector4)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector3,testVector5)):
			selectionFinished = False
		elif (not fourVectors.check_different_direction(testVector4,testVector5)):
			selectionFinished = False
		else:
			selectionFinished = True
	##Have test chain of: q-g-g-g-q_bar.
	possibleParticleNos = [1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,21]
	testParticle1 = particles.particle(possibleParticleNos[random.randrange(0,12)],testVector1)
	testParticle5 = particles.particle(-1*testParticle1.get_code(),testVector5) ##1 not float as code must be an integer.
	##Fill the middle with gluons.
	testParticle2 = particles.particle(21,testVector2)
	testParticle3 = particles.particle(21,testVector3)
	testParticle4 = particles.particle(21,testVector4)
	pC = counters.counter(1)
	testParticle1.set_unique_ID(pC.next())
	testParticle2.set_unique_ID(pC.next())
	testParticle3.set_unique_ID(pC.next())
	testParticle4.set_unique_ID(pC.next())
	testParticle5.set_unique_ID(pC.next())
	testParticles = [testParticle1,testParticle2,testParticle3,testParticle4,testParticle5]
	testRC1 = resultsContainer()
	testRC2 = resultsContainer()
	testRC2.store(testParticle2)
	testRC2.store(testParticle3)
	testRC3 = resultsContainer()
	testCopyRC = testRC2.copy()
	testRCs = [testRC1,testRC2,testRC3]

	##Test check_is_resultsContainer:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_resultsContainer:"
	print "Calling check_is_resultsContainer on instance: " , check_is_resultsContainer(testRC1)
	print "Calling check_is_resultsContainer on wrong instance: " , check_is_resultsContainer(testParticle1)
	print "Calling check_is_resultsContainer on 1.055: " , check_is_resultsContainer(1.055)
	print "Calling check_is_resultsContainer on 'word': " , check_is_resultsContainer('word')
	testResults = [check_is_resultsContainer(testRC1),check_is_resultsContainer(testParticle1)]
	testResults += [check_is_resultsContainer(1.055),check_is_resultsContainer('word')]
	if (sum(testResults) == 1):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_resultsContainer."
	assertions.pause(__name__)

	##Test resultsContainer class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////////////"
	print "Testing resultsContainer class:"
	print "///////////////////////////////"
	assertions.pause(__name__)

	##__init__() tested implicitly.

	##Test store(), get_all() and __getitem__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing store(), get_all() and __getitem__() functions:\n"
	print "The test particles are:"
	for testparticle in testParticles:
		print testparticle.simple_str()
	assertions.pause(__name__)
	print "\nInitially get_all returns:", testRC1.get_all()
	assertions.pause(__name__)
	for testParticle in testParticles:
		testRC1.store(testParticle)
	print "\nAfter storing the test particles it returns:\n", testRC1.get_all()
	testResults = [(testParticle1 == testRC1.get_all()[0]),(testParticle2 == testRC1.get_all()[1])]
	testResults += [(testParticle3 == testRC1.get_all()[2]),(testParticle4 == testRC1.get_all()[3])]
	testResults += [(testParticle5 == testRC1.get_all()[4])]
	if (sum(testResults) == 5):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing store(), get_all() and __getitem__() functions."
	assertions.pause(__name__)

	##Test ___eq__() and __ne__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing __eq__() and __ne__()  functions:\n"
	print "Test results container 1 is:", testRC1.simple_str()
	print "\nTest results container 2 is:", testRC2.simple_str()
	print "\nTest results container 3 is:", testRC3.simple_str()
	print "\n1==1 gives:", testRC1 == testRC1
	print "1==2 gives:", testRC1 == testRC2
	print "1==3 gives:", testRC1 == testRC3
	print "2==2 gives:", testRC2 == testRC2
	print "2==3 gives:", testRC2 == testRC3
	print "3==3 gives:", testRC3 == testRC3
	assertions.pause(__name__)
	print "1!=1 gives:", testRC1 != testRC1
	print "1!=2 gives:", testRC1 != testRC2
	print "1!=3 gives:", testRC1 != testRC3
	print "2!=2 gives:", testRC2 != testRC2
	print "2!=3 gives:", testRC2 != testRC3
	print "3!=3 gives:", testRC3 != testRC3
	testResults = [testRC1 == testRC1,testRC1 == testRC2,testRC1 == testRC3,testRC2 == testRC2,testRC2 == testRC3,testRC3 == testRC3]
	testResults += [testRC1 != testRC1,testRC1 != testRC2,testRC1 != testRC3,testRC2 != testRC2,testRC2 != testRC3,testRC3 != testRC3]
	if (sum(testResults) == 6 and (testRC1 == testRC1)):
		print "\nTest successful!"
	else:
		print "\nTest Failed!"
	print "\nFinished testing __eq__() and __ne__()  functions."
	assertions.pause(__name__)

	##Test __add__() and __iadd__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing __add__() and __iadd__() functions:\n"
	print "Test results container 1 is:\n", testRC1.simple_str()
	print "Test results container 2 is: \n", testRC2.simple_str()
	print "Test results container 3 is: \n", testRC3.simple_str()
	print "Adding 1 and 2 gives:", (testRC1 + testRC2).simple_str()
	print "Adding all three gives:", (testRC1 + testRC2 + testRC3).simple_str()
	print "Adding 1 to itself returns:"
	try:
		print (testRC1 + testRC1).simple_str()
	except(AssertionError):
		print "\nCaused an AssertionError as expected!"
		print "Test successful!"
	print "\nFinished testing __add__() and __iadd__() functions."
	assertions.pause(__name__)

	##Test copy() function:##
	print "\n--------------------------------------------------\n"
	print "Testing copy() function:\n"
	print "Test results container 2 is: \n", testRC2.simple_str()
	print "It's copy is: \n", testCopyRC.simple_str()
	print "Adding them returns:"
	try:
		print testRC2 + testCopyRC
	except(AssertionError):
		print "\nCaused an AssertionError as expected!"
		print "Test successful!"
	print "\nFinished testing copy() function."
	assertions.pause(__name__)

	##Test __str__(), simple_str() and __repr__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing __str__(), simple_str() and __repr__() functions:\n"
	print "Calling in order on the three test results containers used above returns:"
	for testRC in testRCs:
		print "\n", str(testRC)
		print "\n", testRC.simple_str()
		print "\n", repr(testRC)
		assertions.pause(__name__)
	print "\nFinished testing __str__(), simple_str() and __repr__() functions."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "///////////////////////////////////////////"
	print "Finished checking resultsContainers module!"
	print "///////////////////////////////////////////"