####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling dipoles in dipole showering."""

##Import required modules:##
import random
import assertions
import particleData
import particles

print "\n///////////////////////"
print "Loading dipoles module:"
print "///////////////////////\n"

##Functions:##

def check_is_dipole(toCheck):
	"""A function to check for an instance of the dipole class."""
	return isinstance(toCheck,dipole)

##Classes:##

class dipole(object):
	"""A class to handle a dipole."""
	##updated tracks changes in a particle and relevant variables.

	def update_dipole_COM_four_vector(self):
		"""A function to update the COM four-vector of the dipole."""
		__newCOMFourVector = self.__dipoleList[0].get_four_momentum().copy() + self.__dipoleList[1].get_four_momentum().copy()
		self.__COMFourVector = __newCOMFourVector

	def update_dipole_mass_squared(self):
		"""A function to update the mass of the dipole."""
		self.update_dipole_COM_four_vector()
		__COMFourVector = self.__COMFourVector.copy()
		__newDipoleMassSquared = __COMFourVector * __COMFourVector
		self.__massSquared = __newDipoleMassSquared

	def __init__(self,particle1,particle2):
		"""A function to initialise a dipole."""
		##Stored and accessed in the order given to enable order management in chains.
		##Can't have a 'zero dipole' so no check present; creation requires both particles non-zero.
		assert (particles.check_is_particle(particle1) and particles.check_is_particle(particle2))
		assert (particle1.__nonzero__() and particle2.__nonzero__())
		__particle1, __particle2 = particle1.copy(), particle2.copy()
		self.__dipoleList = [__particle1,__particle2]
		self.update_dipole_mass_squared() #Set mass squared and also COM four-vector.
		self.__updated = True

	def set_not_updated(self):
		"""A function to tell a dipole that it needs updating."""
		self.__updated = False

	def copy(self):
		"""A function to create a copy of a dipole."""
		__aCopy = dipole(self.__dipoleList[0].copy(),self.__dipoleList[1].copy())
		if (not self.__updated):
			__aCopy.set_not_updated()
		return __aCopy

	def __getitem__(self,listIndexOfParticle):
		"""A function to return one of the particles in the dipole using its dipole-list index."""
		##No set_particle as update used instead.
		assert assertions.valid_dipole_index(listIndexOfParticle)
		return self.__dipoleList[listIndexOfParticle].copy()

	def __setitem__(self,listIndexToUpdate,updatedParticle):
		"""A function to update a dipole, e.g after a recoil from another dipole has occured."""
		assert assertions.valid_dipole_index(listIndexToUpdate)
		assert particles.check_is_particle(updatedParticle)
		assert updatedParticle.__nonzero__()
		self.__dipoleList[listIndexToUpdate] = updatedParticle.copy()
		self.update_dipole_mass_squared()
		self.__updated = True

	def get_COM_four_vector(self):
		"""A function to return the COM four-vector associated with the dipole."""
		return self.__COMFourVector

	def get_mass_squared(self):
		"""A function to return the mass squared of the dipole."""
		return self.__massSquared

	def get_index_to_recoil(self,processCode):
		"""A function to return the index of the particle in the dipole which would take the recoil during splitting."""
		##Returns 0 or 1
		assert ((1 <= processCode) and (processCode <= 3))
		if (processCode == 1):
			return particles.which_particle_to_recoil_g_emission(self.__dipoleList[0],self.__dipoleList[1])
		elif (processCode == 2):
			__code0, __code1 = self.__dipoleList[0].get_code(), self.__dipoleList[1].get_code()
			__kQs = particleData.knownParticles.get_known_quarks()
			__gC = particleData.knownParticles.get_code_from_name('gluon')
			__expectedQCodes = __kQs + [-x for x in __kQs]
			if ((__code0 in __expectedQCodes) and (__code1 == __gC)):
				return 0
			elif ((__code0 == __gC) and (__code1 in __expectedQCodes)):
				return 1
			elif ((__code0 == __gC) and (__code1 == __gC)):
				__R = random.random() ##Choose randomly as sudakov doesn't specify.
				if (__R < 0.5):
					return 0
				else:
					return 1
		elif (processCode == 3): ##Photon emission only occures between two q(Bar)/q(Bar).
			##For qqBar the cross section is the same (ignoring prefactors) as that of gluon emission.
			##Therefore use the same rules.
			return particles.which_particle_to_recoil_g_emission(self.__dipoleList[0],self.__dipoleList[1])

	def is_updated(self):
		"""A function to see if the dipole is updated."""
		return self.__updated

	def __str__(self):
		"""A function to return a string of a dipole object."""
		__stringOfDipole = ("\n[[ [" + self.__dipoleList[0].get_name() + " - " + self.__dipoleList[1].get_name() + "] | Updated: ")
		__stringOfDipole += (str(self.__updated) + " |\n" + str(self.__dipoleList[0]) + " |\n" + str(self.__dipoleList[1]))
		__stringOfDipole += (" |\n   Mass squared: " + str(self.__massSquared) + "]]")
		return __stringOfDipole

	def simple_str(self):
		"""A function to return a simplified string of a dipole object."""
		__simpleName1 = particleData.get_short_particle_name(self.__dipoleList[0].get_name())
		__simpleName2 = particleData.get_short_particle_name(self.__dipoleList[1].get_name())
		__simpleStringOfDipole = ("[" + __simpleName1 + " - " + __simpleName2 + "]")
		return __simpleStringOfDipole

	def __repr__(self):
		"""A function to return a representation of a dipole object."""
		__reprOfDipole = ("\n<dipole([" + self.__dipoleList[0].get_name() + " - " + self.__dipoleList[1].get_name() + "]| Updated: ")
		__reprOfDipole += (str(self.__updated) + " |\n" + repr(self.__dipoleList[0]) + " |\n" + repr(self.__dipoleList[1]))
		__reprOfDipole += (" |\n   Mass squared: " + str(self.__massSquared) + "])>")
		return __reprOfDipole

	def __eq__(self,other):
		"""A function to check whether two dipoles are equal, up to the code precision."""
		##Just checks if they have identical particles and are updated.
		assert check_is_dipole(other)
		if (self.__dipoleList[0] == other[0]):
			if (self.__dipoleList[1] == other[1]):
				##Require both are updated, as other one might be due a recoil and technically not equal!
				if ((self.__updated == True) and (other.__updated == True)):
					return True
		return False

	def __ne__(self,other):
		"""A function to check whether two dipoles are not equal, up to the code precision."""
		assert check_is_dipole(other)
		return (not self.__eq__(other))

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import math
	import random
	import fourVectors

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////"
	print "Testing dipoles module:"
	print "///////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	##Generate random test four-vectors:##
	testFourVectors1 = {'a':None,'b':None}
	##Prevent the random vectors being the same as can't boost into frame of massless particle.
	testVector1, testVector2 = fourVectors.fourVector(0,0,0,0), fourVectors.fourVector(0,0,0,0) ##Required to start while loop.
	while ((testVector1 == testVector2) or (not fourVectors.check_different_direction(testVector1,testVector2))):
		for testFourVector in testFourVectors1:
			x1 = random.randrange(1,10)/10.0
			x2 = random.randrange(0,int(math.sqrt(1-(x1**2))*10))/10.0
			x3 = random.randrange(0,int(math.sqrt(1-(x2**2))*10))/10.0
			##Make them massless.
			x0 = math.sqrt((x1*x1) + (x2*x2) +(x3*x3))
			testFourVectors1[testFourVector] = fourVectors.fourVector(x0,x1,x2,x3)
		testVector1 = testFourVectors1['a'].copy()
		testVector2 = testFourVectors1['b'].copy()
	##Generate set 2:
	testFourVectors2 = {'a':None,'b':None}
	##Prevent the random vectors being the same as can't boost into frame of massless particle.
	testVector3, testVector4 = fourVectors.fourVector(0,0,0,0), fourVectors.fourVector(0,0,0,0) ##Required to start while loop.
	while ((testVector3 == testVector4) or (not fourVectors.check_different_direction(testVector3,testVector4))):
		for testFourVector in testFourVectors2:
			##Random range set so as to produce a vector slower than the speed of light.
			x1 = random.randrange(1,10)/10.0
			x2 = random.randrange(0,int(math.sqrt(1-(x1**2))*10))/10.0
			x3 = random.randrange(0,int(math.sqrt(1-(x2**2))*10))/10.0
			##Make them massless.
			x0 = math.sqrt((x1*x1) + (x2*x2) +(x3*x3))
			testFourVectors2[testFourVector] = fourVectors.fourVector(x0,x1,x2,x3)
		testVector3 = testFourVectors2['a'].copy()
		testVector4 = testFourVectors2['b'].copy()
	##Test code here doesn't care about matching appropriate particles in the test dipoles.
	possibleParticleNos = [1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,21]
	testParticle1 = particles.particle(possibleParticleNos[random.randrange(0,13)],testVector1)
	testParticle2 = particles.particle(possibleParticleNos[random.randrange(0,13)],testVector2)
	testParticle3 = particles.particle(possibleParticleNos[random.randrange(0,13)],testVector3)
	testParticle4 = particles.particle(possibleParticleNos[random.randrange(0,13)],testVector4)
	testDipole1 = dipole(testParticle1,testParticle2)
	testDipole2 = dipole(testParticle3,testParticle4)
	testCopyDipole2 = testDipole2.copy()
	testCopy2Dipole1 = testDipole1.copy()
	testCopy2Dipole1.set_not_updated()
	listOfCodes = [1,2,3,4,5,6,-1,-2,-3,-4,-5,-6,21]

	##Test check_is_dipole:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_dipole:\n"
	print "Calling check_is_dipole on instance: " , check_is_dipole(testDipole1)
	print "Calling check_is_dipole on wrong instance: " , check_is_dipole(testVector1)
	print "Calling check_is_dipole on 1.055: " , check_is_dipole(1.055)
	print "Calling check_is_dipole on 'word': " , check_is_dipole('word')
	results = [check_is_dipole(testDipole1),check_is_dipole(testVector1),check_is_dipole(1.055),check_is_dipole('word')]
	if (not ((sum(results) != 1) or (results[0] != True))):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_dipole."
	assertions.pause(__name__)

	##Test dipole class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////"
	print "Testing dipole class:"
	print "/////////////////////"
	assertions.pause(__name__)

	##__init__ tested implicitly.

	##Test __str__(), simple_str() and __repr__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing __str__(), simple_str() and __repr__() functions:\n"
	print "Starting with test dipole1:\n", str(testDipole1)
	assertions.pause(__name__)
	print "\nThe simplified string gives:\n", testDipole1.simple_str()
	print "\n__repr__() gives:\n", repr(testDipole1)
	print "\nFinished testing __str__(), simple_str() and __repr__() functions."
	assertions.pause(__name__)

	##Test get_COM_four_vector(), get_mass_squared(), set_not_updated() and is_updated().##
	print "\n--------------------------------------------------\n"
	print "Testing get_COM_four_vector(), get_mass_squared(), set_not_updated() and is_updated().\n"
	print "Looking at testDipole1:\n" , testDipole1
	print "\nget_COM_four_vector return:", testDipole1.get_COM_four_vector()
	print "and get_mass_squared returns:", testDipole1.get_mass_squared()
	print "Calling is_updated returns:", testDipole1.is_updated()
	testDipole1.set_not_updated()
	print "Calling set_not_updated and then is_updated returns:", testDipole1.is_updated()
	assertions.pause(__name__)
	print "\nLooking at testDipole2:\n" , testDipole2
	print "\nget_COM_four_vector return:", testDipole2.get_COM_four_vector()
	print "and get_mass_squared returns:", testDipole2.get_mass_squared()
	print "Calling is_updated returns:", testDipole2.is_updated()
	testDipole2.set_not_updated()
	print "Calling set_not_updated and then is_updated returns:", testDipole2.is_updated()
	if not (testDipole1.is_updated() and testDipole2.is_updated()):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_COM_four_vector(), get_mass_squared(), set_not_updated() and is_updated()."
	assertions.pause(__name__)

	##Test  copy(), __getitem__() and __setitem__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing  copy(), __getitem__() and __setitem__() functions:\n"
	print "Starting with test dipole:\n", testDipole1
	testCopyDipole1 = testDipole1.copy()
	print "\nCreating a copy gives:\n", testCopyDipole1
	assertions.pause(__name__)
	print "\nCalling __getitem__ on both elements of dipole1 gives:\n", testCopyDipole1[0], "\n", testCopyDipole1[1]
	print "and on element 1 of dipole2 gives:\n", testCopyDipole2[0]
	assertions.pause(__name__)
	testCopyDipole1[0] = testCopyDipole2[0]
	print "\nThen setting particle1 of dipole1 to that of dipole2 gives:\n", testCopyDipole1
	if (testCopyDipole1[0] == testCopyDipole2[0]):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing  copy(), __getitem__() and __setitem__() functions."
	assertions.pause(__name__)

	##Test  update_dipole_COM_four_vector and update_dipole_mass_squared functions:##
	print "\n--------------------------------------------------\n"
	print "Testing  update_dipole_COM_four_vector and update_dipole_mass_squared functions:\n"
	print "These functions should automatically be run when a particle is updated via __setitem__."
	print "For dipole1 we had:\n", testDipole1.get_COM_four_vector(), "\nand:", testDipole1.get_mass_squared()
	print "Now having updated particle1 of the copy it has:\n", testCopyDipole1.get_COM_four_vector(), "\nand:", testCopyDipole1.get_mass_squared()
	test1 = (testCopyDipole1.get_COM_four_vector() != testDipole1.get_COM_four_vector())
	test2 = (testCopyDipole1.get_mass_squared() != testDipole1.get_mass_squared())
	if (test1 and test2):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing  update_dipole_COM_four_vector and update_dipole_mass_squared functions."
	assertions.pause(__name__)

	##Test __eq__() and __ne__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing __eq__() and __ne__() functions:\n"
	print "Using testDipole1:\n", testDipole1
	print "\nand testCopyDipole1:\n", testCopyDipole1
	assertions.pause(__name__)
	print "Two dipoles can't be equal if updated == False as for testDipole1:"
	print "testDipole1 == testDipole1:", (testDipole1 == testDipole1)
	testDipole1._dipole__updated = True ##WARNING: Name mangling!###
	print "\nNow setting it to True:\n"
	print "testDipole1 == testDipole1:", (testDipole1 == testDipole1)
	print "testCopyDipole1 == testCopyDipole1:", (testCopyDipole1 == testCopyDipole1)
	print "testCopyDipole1 == testDipole1:", (testCopyDipole1 == testDipole1)
	print "testDipole1 == testCopyDipole1:", (testDipole1 == testCopyDipole1)
	results = [(testDipole1 == testDipole1),(testCopyDipole1 == testCopyDipole1),(testCopyDipole1 == testDipole1),(testDipole1 == testCopyDipole1)]
	print "testDipole1 != testDipole1:", (testDipole1 != testDipole1)
	print "testCopyDipole1 != testCopyDipole1:", (testCopyDipole1 != testCopyDipole1)
	print "testCopyDipole1 != testDipole1:", (testCopyDipole1!= testDipole1)
	print "testDipole1 != testCopyDipole1:", (testDipole1 != testCopyDipole1)
	results += [(testDipole1!= testDipole1),(testCopyDipole1 != testCopyDipole1),(testCopyDipole1 != testDipole1),(testDipole1 != testCopyDipole1)]
	print "\nTesting against testCopy2Dipole1, == testDipole1 but with updated = False:"
	print "\ntestDipole1 == testCopy2Dipole1:", (testDipole1 == testCopy2Dipole1)
	print "testDipole1 != testCopy2Dipole1:", (testDipole1 != testCopy2Dipole1)
	results += [(testDipole1 == testCopy2Dipole1),(testDipole1 != testCopy2Dipole1)]
	if (not ((sum(results) != 5) or (results[0] != True) or (results[1] != True) or (results[6] != True) or (results[7] != True) or (results[9] != True))):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __eq__() and __ne__() functions."
	assertions.pause(__name__)

	##Test get_index_to_recoil function:##
	##Already tested extensively in particles module.
	print "\n--------------------------------------------------\n"
	print "Testing get_index_to_recoil function:\n"
	for code1 in listOfCodes:
		for code2 in listOfCodes:
			tPN1 = particles.particle(code1,fourVectors.fourVector(2.0,0.0,0.0,0.0))
			tPN2 = particles.particle(code2,fourVectors.fourVector(1.0,0.0,0.0,0.0))
			tDN = dipole(tPN1,tPN2)
			print "Given the dipole:\n", tDN
			print "\n returns:", tDN.get_index_to_recoil()
			assertions.pause(__name__)
	print "\nFinished testing get_index_to_recoil function."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "/////////////////////////////////"
	print "Finished checking dipoles module!"
	print "/////////////////////////////////"