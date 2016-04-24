####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module containing a class to handle a particle throughout dipole showering."""

##Import required modules:##
import random
import assertions
import precision
import particleData
import fourVectors

print "\n/////////////////////////"
print "Loading particles module:"
print "/////////////////////////\n"

##Functions:##

def check_is_particle(toCheck):
	"""A function to check for an instance of the particle class."""
	return isinstance(toCheck,particle)

def which_particle_most_energic(particle1,particle2):
	"""A function to determine which of two particles is most energetic."""
	assert (check_is_particle(particle1) and check_is_particle(particle2))
	assert particle1.__nonzero__() and particle2.__nonzero__()
	##Returns relevant dipole list index.
	if (particle1.get_four_momentum()[0] > particle2.get_four_momentum()[0]):
		return 0
	elif (particle1.get_four_momentum()[0] < particle2.get_four_momentum()[0]):
		return 1
	else: ##If they have the same energy, randomly choose one of them to decay:
		return random.randrange(0,2)

def which_particle_to_recoil_g_emission(particle1,particle2):
	"""A function to determine which of the two particles from a dipole will take the recoil if they emit a gluon."""
	assert (check_is_particle(particle1) and check_is_particle(particle2))
	##Source: Paper 'Ariadne version 4 - A program for simulation of QCD cascades implementing the colour dipole model'.
	##Returns relevant dipole list index.
	##calc_f_for_qg in sudakov is written (tX3^3 not tX1^3) under these assumptions.
	##q-g, g-q, qBar-g or g-qBar -> g recoils:
	if (particle1.is_quark() and particle2.get_name() == "gluon"):
		return 1
	elif (particle1.get_name() == "gluon" and particle2.is_quark()):
		return 0
	##g-g -> Most energetic recoils:
	elif ((particle1.get_name() == "gluon") and (particle1.get_name() == "gluon")):
		#Ariadne does something angle related here?
		return which_particle_most_energic(particle1,particle2)
	##Source: 'From two to three jets in heavy boson decays: An algorithmic approach'.
	##q-q_bar -> Let one quark retain its direction with a probability proportional to its energy^2:
	elif (particle1.is_quark() and particle2.is_quark()):
		__E1, __E2 = particle1.get_four_momentum()[0], particle2.get_four_momentum()[0]
		__E1Squared, __E2Squared = __E1*__E1, __E2*__E2
		__summedEnergySquared = __E1Squared + __E2Squared
		if (random.uniform(0.0,__summedEnergySquared) <= __E1Squared):
			##Let quark 1 retain its direction:
			return 1
		else:
			##Let quark 2 retain its direction:
			return 0

##Classes:##

class particle(object):
	"""A class for a particle."""

	def __init__(self,lookUpCodeOrName,fourMomentum = fourVectors.fourVector(),mothers = [0,0],children = [0,0],colours = [0,0], statusCode = 1):
		"""A function to initiate a particle."""
		##Allow initiation via name or code but then works from the code.
		##StatusCode: -1 initial, 0 = intermediate, 1 = final state. 1 chosen by default as this is currently a final state shower.
		assert ((type(lookUpCodeOrName) == int) or (type(lookUpCodeOrName) == str))
		assert fourVectors.check_is_fourVector(fourMomentum)
		assert ((type(mothers) == list) and (type(mothers[0]) == int) and (type(mothers[1]) == int))
		assert ((type(children) == list) and (type(children[0]) == int) and (type(children[1]) == int))
		assert ((type(colours) == list) and (type(colours[0]) == int) and (type(colours[1]) == int))
		assert (type(statusCode) == int)
		assert assertions.valid_particle_status_code(statusCode)
		if (type(lookUpCodeOrName) == int):
			assert particleData.knownParticles.has_key(lookUpCodeOrName)
			self.__lookUpCode = lookUpCodeOrName
		elif (type(lookUpCodeOrName) == str):
			assert particleData.knownParticles.has_name(lookUpCodeOrName)
			self.__lookUpCode = particleData.knownParticles.get_code_from_name(lookUpCodeOrName)
		self.__fourMomentum = fourMomentum
		self.__mothers = mothers
		self.__children = children
		self.__colours = colours
		self.__statusCode = statusCode
		self.__historyList = [False,False,False,False] ##Mothers then children (if set)
		self.__colourSetList = [False, False]
		for __i, __a in enumerate(self.__mothers+self.__children): ##Do both at once.
			if __a != 0:
				self.__historyList[__i] = True
		for __i, __a in enumerate(self.__colours):
			if __a != 0:
				self.__colourSetList[__i]= True
		self.__uniqueID = 0
		self.__producedAt = -1 ##Default value for not set, as 0 is actually used.

	def __nonzero__(self):
		"""A function to determine whether the particle's four-momentum is a zero object (empty / not fully initiated)."""
		return self.__fourMomentum.__nonzero__()

	def __getitem__(self,index):
		"""A function to get a component of the energy-momentum fourVector of a particle."""
		assert assertions.all_are_numbers([index])
		assert assertions.valid_four_vector_index(index)
		return self.__fourMomentum[index]

	def get_four_momentum(self):
		"""A function to return the four-momentum of a particle."""
		return self.__fourMomentum.copy()
      
	def set_four_momentum(self,fourMomentumToSet):
		"""A function to set the four-momentum of a particle."""
		assert fourVectors.check_is_fourVector(fourMomentumToSet)
		assert fourMomentumToSet.__nonzero__()
		self.__fourMomentum = fourMomentumToSet.copy()

	def get_unique_ID(self):
		"""A function to get a particle's unique ID."""
		return self.__uniqueID

	def set_unique_ID(self,toSet,overwrite = 0):
		"""A function to set a particles unique ID."""
		assert (type(toSet) == int)
		assert ((type(overwrite) == int) and ((overwrite == 0) or (overwrite == 1)))
		if (overwrite == 1): ##To prevent errors: Must specify overwriting if already set.
			self.__uniqueID = toSet
		else:
			if (self.__uniqueID == 0):
				self.__uniqueID = toSet
			else:
				##Don't overwrite if already set and not specified to.
				print "Trying to overwrite ID!"
				assert False

	def get_produced_at(self):
		"""A function to get the value of Pperp squared a particle was produced at."""
		return self.__producedAt

	def set_produced_at(self,toSet,overwrite = 0):
		"""A function to set the value of Pperp squared a particle was produced at."""
		assert (type(toSet) == float or (toSet == -1))
		assert ((type(overwrite) == int) and ((overwrite == 0) or (overwrite == 1)))
		if (overwrite == 1): ##To prevent errors: Must specify overwriting if already set.
			self.__producedAt = toSet
		else:
			if (self.__producedAt == -1):
				self.__producedAt = toSet
			else:
				##Don't overwrite if already set and not specified to.
				print "Trying to overwrite producedAt!"
				assert False

	def copy(self):
		"""A function to create an exact copy of a particle."""
		__new = particle(self.__lookUpCode,self.__fourMomentum.copy(),self.__mothers,self.__children,self.__colours,self.__statusCode)
		__new.set_unique_ID(self.__uniqueID)
		__new.set_produced_at(self.__producedAt)
		return __new

	def __str__(self):
		"""A function to return a string of a particle object."""
		__theStr = ("[[" + self.get_name() + " | " + str(self.__fourMomentum) + " |\n" + str(self.__mothers) + " " + str(self.__children))
		__theStr += (" | " + str(self.__colours) + " | " + str(self.__statusCode))
		__theStr += (" | " + str(self.__uniqueID) + " " + str(self.__producedAt) + " ]]")
		return __theStr

	def simple_str(self):
		"""A function to return a simple string of a particle."""
		return "[[" + str(self.__uniqueID) + " : " + str(self.get_name()) + "]]"

	def __repr__(self):
		"""A function to return a representation of a particle object."""
		__theRepr = ("<particle([" + self.get_name() + " |\n" + repr(self.__fourMomentum) + " |\n" + str(self.__mothers) + " " + str(self.__children))
		__theRepr += (" | " + str(self.__colours) + " | " + str(self.__statusCode))
		__theRepr += (" | " + str(self.__uniqueID) + " " + str(self.__producedAt) + " ])>")
		return __theRepr

	def __eq__(self,other):
		"""A function to check whether two particles are equal up to the code precision."""
		assert check_is_particle(other)
		assert (self.__nonzero__() and other.__nonzero__())
		if (self.__lookUpCode == other.get_code()):
			if (self.__fourMomentum == other.get_four_momentum()):
				__checks = [(self.get_mother(0) == other.get_mother(0)),(self.get_mother(1) == other.get_mother(1))]
				__checks += [(self.get_child(0) == other.get_child(0)),(self.get_child(1) == other.get_child(1))]
				__checks += [(self.get_colour(0) == other.get_colour(0)),(self.get_colour(1) == other.get_colour(1))]
				__checks += [(self.get_status_code() == other.get_status_code()),(self.get_unique_ID() == other.get_unique_ID())]
				__checks += [(self.get_produced_at() == other.get_produced_at())]
				if (sum(__checks) == 9):
					return True
		return False

	def __ne__(self,other):
		"""A function to check whether two particles are not equal up to the code precision."""
		return not self.__eq__(other)

	def get_code(self):
		"""A function to return the PDG code of the particle."""
		return self.__lookUpCode

	def get_name(self):
		"""A function to return the name of the particle."""
		return particleData.knownParticles.get_name_from_code(self.__lookUpCode)

	def get_mass_from_code(self):
		"""A function to return the normal rest mass of the particle."""
		##Uses the PDG value. Not correct when working with massless particles.
		return particleData.knownParticles.get_mass_from_code(self.__lookUpCode)

	def get_mass(self):
		"""A function to return the mass of the particle from its four-momentum."""
		##Mass is Minkowski scalar product of the E-p four-vector.
		assert self.__nonzero__()
		__v = self.get_four_momentum()
		return __v*__v

	def get_width(self):
		"""A function to return the width of the particle."""
		return particleData.knownParticles.get_width_from_code(self.__lookUpCode)

	def get_charge(self):
		"""A function to return the charge of the particle."""
		return particleData.knownParticles.get_charge_from_code(self.__lookUpCode)

	def get_spin(self):
		"""A function to return the spin of the particle."""
		return particleData.knownParticles.get_spin_from_code(self.__lookUpCode)

	def get_mother(self,index):
		"""A function to return a mother of the particle given its index."""
		assert (type(index) == int)
		assert assertions.valid_dipole_index(index) ##As a list of the same dimensions.
		return self.__mothers[index]

	def set_mother(self,index,newValue):
		"""A function to update a mother of the particle given its index and a new value."""
		assert (type(index) == int)
		assert assertions.valid_dipole_index(index) ##As a list of the same dimensions.
		assert (type(newValue) == int)
		self.__mothers[index] = newValue
		if (newValue != 0):
			self.__historyList[index] = True
		else:
			self.__historyList[index] = False

	def get_child(self,index):
		"""A function to return a child of the particle given its index."""
		assert (type(index) == int)
		assert assertions.valid_dipole_index(index) ##As a list of the same dimensions.
		return self.__children[index]

	def set_child(self,index,newValue):
		"""A function to update a child of the particle given its index and a new value."""
		assert (type(index) == int)
		assert assertions.valid_dipole_index(index) ##As a list of the same dimensions.
		assert (type(newValue) == int)
		self.__children[index] = newValue
		if (newValue != 0):
			self.__historyList[index+2] = True
		else:
			self.__historyList[index+2] = False

	def get_colour(self,index):
		"""A function to return a colour of the particle given its index."""
		assert (type(index) == int)
		assert assertions.valid_dipole_index(index) ##As a list of the same dimensions.
		return self.__colours[index]

	def set_colour(self,index,newValue):
		"""A function to update a colour of the particle given its index and a new value."""
		assert (type(index) == int)
		assert assertions.valid_dipole_index(index) ##As a list of the same dimensions.
		assert (type(newValue) == int)
		self.__colours[index] = newValue
		if newValue != 0:
			self.__colourSetList[index] = True
		else:
			self.__colourSetList[index] = False

	def get_status_code(self):
		"""A function to return the status code of the particle."""
		return self.__statusCode

	def is_lepton(self):
		"""A function to check if the particle is a lepton."""
	 	return particleData.knownParticles.is_lepton(self.__lookUpCode)

	def is_quark(self):
		"""A function to check if the particle is a quark."""
		return particleData.knownParticles.is_quark(self.__lookUpCode)

	def is_meson(self):
		"""A function to check if the particle is a meson."""
		return particleData.knownParticles.is_meson(self.__lookUpCode)

	def is_hadron(self):
		"""A function to check if the particle is a hadron."""
		return particleData.knownParticles.is_hadron(self.__lookUpCode)

	def is_baryon(self):
		"""A function to check if the particle is a baryon."""
		return particleData.knownParticles.is_baryon(self.__lookUpCode)

	def is_fermion(self):
		"""A function to check if the particle is a fermion."""
		return particleData.knownParticles.is_fermion(self.__lookUpCode)

	def is_gauge_boson(self):
		"""A function to check if the particle is a gauge boson."""
		return particleData.knownParticles.is_gauge_boson(self.__lookUpCode)

	def is_boson(self):
		"""A function to check if the particle is a boson."""
		return particleData.knownParticles.is_boson(self.__lookUpCode)

	def is_anti_particle(self):
		"""A function to check if the particle is an anti-particle."""
		return particleData.knownParticles.is_anti_particle(self.__lookUpCode)

	def history_set(self):
		"""A function to return whether the history of a particle is set."""
		if (sum(self.__historyList) == 4):
			return True
		else:
			return False

	def colour_set(self):
		"""A function to return whether a colour for a particle is set."""
		if (sum(self.__colourSetList) >= 1): ##Currently only checks that one colour is set.
			return True
		else:
			return False

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	from matplotlib import pyplot
	from matplotlib import lines

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////"
	print "Testing particles module:"
	print "////////////////////////"
	assertions.pause(__name__)

	##Setup here:##
	print "\nGenerating test values..."
	##Generate random test four-vectors:##
	tFourVectors = {'a':None,'b':None,'c':None,'d':None}
	for tFourVector in tFourVectors:
		tX0 = random.randrange(20)
		tX1 = random.randrange(20)
		tX2 = random.randrange(20)
		tX3 = random.randrange(20)
		tFourVectors[tFourVector] = fourVectors.fourVector(tX0,tX1,tX2,tX3)
	tZeroParticle = particle(3)
	tParticle1 = particle(2,tFourVectors['a'].copy(),[5,6],[7,8],[51,52],-1)
	tParticle2 = particle(1,tFourVectors['b'].copy(),[0,0],[0,0],[0,0],-1)
	tParticle3 = particle('photon',tFourVectors['c'].copy())
	tParticle4 = tZeroParticle.copy()
	tParticle4.set_four_momentum(tParticle1.get_four_momentum())
	tP1C = tParticle1.copy()
	tParticle5 = particle(2,tFourVectors['a'].copy(),[5,6],[7,8],[51,52],-1)
	tQuark1 = particle(1,tFourVectors['a'].copy())
	tQuark2FourVector = tFourVectors['a'].copy()
	tQuark2FourVector *= 2
	tGluon2FourVector = tQuark2FourVector.copy()
	tQuark2 = particle(2,tQuark2FourVector)
	tGluon1 = particle(21,tFourVectors['a'].copy())
	tGluon2 = particle(21,tGluon2FourVector)
	tNumTests = 1000

	##Test check_is_particle:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_particle:\n"
	tSuccessful = True
	print "Calling on tParticle1:" , check_is_particle(tParticle1)
	if not check_is_particle(tParticle1):
		tSuccessful = False
	print "Calling on tFourVectors['a']:" , check_is_particle(tFourVectors['a'])
	if check_is_particle(tFourVectors['a']):
		tSuccessful = False
	print "Calling on 1.055:" , check_is_particle(1.055)
	if check_is_particle(1.055):
		tSuccessful = False
	print "Calling on 'word':" , check_is_particle('word')
	if check_is_particle('word'):
		tSuccessful = False
	print "Calling on tParticle2:" , check_is_particle(tParticle2)
	if not check_is_particle(tParticle2):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_particle."
	assertions.pause(__name__)

	##Test which_particle_most_energic:##
	print "\n--------------------------------------------------\n"
	print "Testing which_particle_most_energic:\n"
	tSuccessful = True
	print "Calling on:", tParticle1, "\n", tParticle2
	print "returns:", which_particle_most_energic(tParticle1,tParticle2)
	tExpectRandom = False
	if tParticle1[0] > tParticle2[0]:
		tExpected1,tExpected2 = 0, 1
	elif tParticle1[0] == tParticle2[0]:
		tExpectRandom = True ##To handle the fact that a random result is returned if the energies match.
	else:
		tExpected1, tExpected2 = 1, 0
	if (not tExpectRandom):
		if (not (which_particle_most_energic(tParticle1,tParticle2) == tExpected1)):
			tSuccessful = False
	print "\nCalling on:", tParticle2, "\n", tParticle1
	print "returns:", which_particle_most_energic(tParticle2,tParticle1)
	if (not tExpectRandom):
		if (not (which_particle_most_energic(tParticle2,tParticle1) == tExpected2)):
			tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	assertions.pause(__name__)
	print "\nTest when the same energy; Should randomly return 0 or 1:"
	print "Calling on:", tParticle1, "\n", tParticle1
	print "returns:", which_particle_most_energic(tParticle1,tParticle1)
	print "\nFinished testing which_particle_most_energic."
	assertions.pause(__name__)

	##Test which_particle_to_recoil_g_emission:##
	print "\n--------------------------------------------------\n"
	print "Testing which_particle_to_recoil_g_emission:\n"
	tSuccessful = True
	print "The q-q_bar (or q-q here) pairs should recoil with a probability weighted by their energy squared."
	print "Calling on:", tQuark1, "\n", tQuark2
	print "returns:", which_particle_to_recoil_g_emission(tQuark1,tQuark2)
	if (which_particle_to_recoil_g_emission(tQuark1,tQuark2) == which_particle_most_energic(tQuark1,tQuark2)):
		print "\nThe more energetic particle chosen to recoil."
	else:
		print "\nThe less energetic particle chosen to recoil."
	print "\nCalling on:", tQuark2, "\n", tQuark1
	print "returns:", which_particle_to_recoil_g_emission(tQuark2,tQuark1)
	if (which_particle_to_recoil_g_emission(tQuark2,tQuark1) == which_particle_most_energic(tQuark2,tQuark1)):
		print "\nThe more energetic particle chosen to recoil."
	else:
		print "\nThe less energetic particle chosen to recoil."
	print "\nCalling on:", tGluon1, "\n", tGluon2
	print "returns:", which_particle_to_recoil_g_emission(tGluon1,tGluon2)
	if (not (which_particle_to_recoil_g_emission(tGluon1,tGluon2) == which_particle_most_energic(tGluon1,tGluon2))):
		tSuccessful = False
	assertions.pause(__name__)
	print "\nCalling on:", tGluon2, "\n", tGluon1
	print "returns:", which_particle_to_recoil_g_emission(tGluon2,tGluon1)
	if (not (which_particle_to_recoil_g_emission(tGluon2,tGluon1) == which_particle_most_energic(tGluon2,tGluon1))):
		tSuccessful = False
	print "\nCalling on:", tQuark1, "\n", tGluon1
	print "returns:", which_particle_to_recoil_g_emission(tQuark1,tGluon2)
	if (not (which_particle_to_recoil_g_emission(tQuark1,tGluon2) == 1)):
		tSuccessful = False
	print "\nCalling on:", tGluon2, "\n", tQuark2
	print "returns:", which_particle_to_recoil_g_emission(tGluon1,tQuark2)
	if (not (which_particle_to_recoil_g_emission(tGluon1,tQuark2) == 0)):
		tSuccessful = False
	if tSuccessful:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	assertions.pause(__name__)
	print "\nTest each particle with itself should give randomly chosen results:"
	print "\nCalling on:", tQuark1, "\n", tQuark1
	print "returns:", which_particle_to_recoil_g_emission(tQuark1,tQuark1)
	print "\nCalling on:", tQuark2, "\n", tQuark2
	print "returns:", which_particle_to_recoil_g_emission(tQuark2,tQuark2)
	print "\nCalling on:", tGluon1, "\n", tGluon1
	print "returns:", which_particle_to_recoil_g_emission(tGluon1,tGluon1)
	print "\nCalling on:", tGluon2, "\n", tGluon2
	print "returns:", which_particle_to_recoil_g_emission(tGluon2,tGluon2)
	print "\nFinished testing which_particle_to_recoil_g_emission."
	assertions.pause(__name__)

	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////"
	print "Testing particle class:"
	print "////////////////////////"
	assertions.pause(__name__)

	##__init__ tested implicitly
	##__getitem__ tested implicitly when testing the functions above.

	##Test __nonzero__, copy, __str__, simple_str, get_four_momentum and set_four_momentum:##
	print "\n--------------------------------------------------\n"
	print "Testing __nonzero__, copy, __str__, simple_str, get_four_momentum and set_four_momentum:\n"
	print "Calling __nonZero__ on tZeroParticle: " , tZeroParticle.__nonzero__()
	print "Calling __nonZero__ on tParticle1: " , tParticle1.__nonzero__()
	print "The zero particle is:\n" , tZeroParticle
	print "\ntParticle1 is:\n" , tParticle1
	print "\nIts simple string is:", tParticle1.simple_str()
	print "\nCoying tZeroParticle and setting four-momentum to tParticle1's gives:"
	print tParticle4
	print "\nWhilst tZeroParticle is still:\n" , tZeroParticle
	tP1C = tParticle1.copy()
	print "\nCoying tParticle1 gives:"
	print tP1C
	tResults = [tZeroParticle.__nonzero__(),tParticle1.__nonzero__(),tParticle4.__nonzero__()]
	tResults += [(tP1C.get_mother(0) == tParticle1.get_mother(0)),(tP1C.get_mother(1) == tParticle1.get_mother(1))]
	tResults += [(tP1C.get_child(0) == tParticle1.get_child(0)),(tP1C.get_child(1) == tParticle1.get_child(1))]
	tResults += [(tP1C.get_colour(0) == tParticle1.get_colour(0)),(tP1C.get_colour(1) == tParticle1.get_colour(1))]
	tResults += [(tP1C.get_status_code() == tParticle1.get_status_code())]
	if (sum(tResults) == 9):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __nonzero__, copy, __str__, simple_str, get_four_momentum and set_four_momentum."
	assertions.pause(__name__)

	##Test get_, set_unique_ID and get_, set_produced_at:##
	tResults = []
	print "\n--------------------------------------------------\n"
	print "Testing get_, set_unique_ID and get_, set_produced_at:\n"
	print "Using:\n", tParticle5
	print "\nCalling get_unique_ID returns:", tParticle5.get_unique_ID()
	print "Calling get_produced_at returns:", tParticle5.get_produced_at()
	tResults += [tParticle5.get_unique_ID(),tParticle5.get_produced_at()]
	print "\nThen call set with 28 and 67.659:"
	tParticle5.set_unique_ID(28)
	tParticle5.set_produced_at(67.659)
	print "\nNow calling get_unique_ID returns:", tParticle5.get_unique_ID()
	print "and calling get_produced_at returns:", tParticle5.get_produced_at()
	tResults += [tParticle5.get_unique_ID(),tParticle5.get_produced_at()]
	print "\nThen try calling set with 113 and 42.555 but without overwrite:"
	try:
		tParticle5.set_unique_ID(113)
		tResults += [False]
	except AssertionError as e:
		print "Trying set_unique_ID produced the error: AssertionError"
		tResults += [True]
	try:
		tParticle5.set_produced_at(42.555)
		tResults += [False]
	except AssertionError as e:
		print "Trying set_produced_at produced the error: AssertionError"
		tResults += [True]
	print "\nThen call set with 113 and 42.555 with overwrite:"
	tParticle5.set_unique_ID(113,1)
	tParticle5.set_produced_at(42.555,1)
	print "\nNow calling get_unique_ID returns:", tParticle5.get_unique_ID()
	print "and calling get_produced_at returns:", tParticle5.get_produced_at()
	tResults += [tParticle5.get_unique_ID(),tParticle5.get_produced_at()]
	tExpected = [0,-1,28,67.659,True,True,113,42.555]
	if (tResults == tExpected):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_, set_unique_ID and get_, set_produced_at."
	assertions.pause(__name__)

	##Test repr, __eq__ and __ne__:##
	print "\n--------------------------------------------------\n"
	print "Testing repr, __eq__ and __ne__:\n"
	print "repr(tZeroParticle) is:\n" , repr(tZeroParticle)
	print "repr(tParticle3) is:\n" , repr(tParticle3)
	print "\n2 == 3 gives:" , tParticle2 == tParticle3
	print "2 == 2 gives:" , tParticle2 == tParticle2
	print "3 == 2 gives:" , tParticle3 == tParticle2
	print "3 == 3 gives:" , tParticle3 == tParticle3
	print "2 != 3 gives:" , tParticle2 != tParticle3
	print "2 != 2 gives:" , tParticle2 != tParticle2
	print "3 != 2 gives:" , tParticle3 != tParticle2
	print "3 != 3 gives:" , tParticle3 != tParticle3
	tResults = [tParticle2 == tParticle3,tParticle2 == tParticle2,tParticle3 == tParticle2]
	tResults += [tParticle3 == tParticle3,tParticle2 != tParticle3,tParticle2 != tParticle2]
	tResults += [tParticle3 != tParticle2,tParticle3 != tParticle3]
	if (sum(tResults) == 4):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing repr, __eq__ and __ne__."
	assertions.pause(__name__)

	##Test the other get__'' functions:##
	print "\n--------------------------------------------------\n"
	print "Testing the other get__'' functions:\n"
	print "Using tParticle2:\n" , tParticle2
	print "get_code returns:" , tParticle2.get_code()
	print "get_name returns:" , tParticle2.get_name()
	print "get_mass_from_code returns:" , tParticle2.get_mass_from_code()
	print "get_mass returns:" , tParticle2.get_mass()
	print "get_width returns:" , tParticle2.get_width()
	print "get_charge returns:" , tParticle2.get_charge()
	print "get_spin returns:", tParticle2.get_spin()
	print "get_mother(0) returns", tParticle2.get_mother(0)
	print "get_mother(1) returns", tParticle2.get_mother(1)
	print "get_child(0) returns", tParticle2.get_child(0)
	print "get_child(1) returns", tParticle2.get_child(1)
	print "get_colour(0) returns", tParticle2.get_colour(0)
	print "get_colour(1) returns", tParticle2.get_colour(1)
	print "get_status_code returns", tParticle2.get_status_code()
	tResults = [tParticle2.get_code(),tParticle2.get_name(),tParticle2.get_mass_from_code(),tParticle2.get_mass()]
	tResults += [tParticle2.get_width(),tParticle2.get_charge(),tParticle2.get_spin(),tParticle2.get_mother(0)]
	tResults += [tParticle2.get_mother(1),tParticle2.get_child(0),tParticle2.get_child(1),tParticle2.get_colour(0)]
	tResults += [tParticle2.get_colour(1),tParticle2.get_status_code()]
	tP2FM = tParticle2.get_four_momentum()
	if (tResults== [1,'d-quark',0.0048,tP2FM*tP2FM,0.0,-1.0/3.0,0.5,0,0,0,0,0,0,-1]):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing the other get__'' functions."
	assertions.pause(__name__)

	##Test colour_set, history_set, set: _mother, _child, _colour'' functions:##
	print "\n--------------------------------------------------\n"
	print "Testing colour_set, history_set, set: _mother, _child, _colour'' functions:\n"
	print "Using tParticle2:\n" , tParticle2
	print "Calling history_set returns:", tParticle2.history_set()
	print "get_mother(0) returns", tParticle2.get_mother(0)
	print "get_mother(1) returns", tParticle2.get_mother(1)
	tParticle2.set_mother(0,9)
	tParticle2.set_mother(1,8)
	print "Setting to 9, 8 and calling again returns:", tParticle2.get_mother(0), tParticle2.get_mother(1)
	print "get_child(0) returns", tParticle2.get_child(0)
	print "get_child(1) returns", tParticle2.get_child(1)
	tParticle2.set_child(0,21)
	tParticle2.set_child(1,22)
	print "Setting to 21, 22 and calling again returns:", tParticle2.get_child(0), tParticle2.get_child(1)
	print "Calling history_set now returns:", tParticle2.history_set()
	print "Calling colour_set returns:", tParticle2.colour_set()
	print "get_colour(0) returns:", tParticle2.get_colour(0)
	print "get_colour(1) returns:", tParticle2.get_colour(1)
	tParticle2.set_colour(0,501)
	tParticle2.set_colour(1,502)
	print "Setting to 501 and 502 and calling again returns:", tParticle2.get_colour(0), tParticle2.get_colour(1)
	print "Calling colour_set now returns:", tParticle2.colour_set()
	if (tParticle2.history_set() and tParticle2.colour_set()):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing colour_set, history_set, set: _mother, _child, _colour'' functions."
	assertions.pause(__name__)

	##Test is__'' functions on all possible particles:##
	print "\n--------------------------------------------------\n"
	print "Testing is__'' functions on all possible particles:\n"
	for tLookUpCode in particleData.knownParticles.get_all_codes():
		if random.random() < 0.5:
			tLookUpCode *= -1 ##Codes must be an integer.
		tParticle = particle(tLookUpCode)
		print "Using tParticle:\n" , tParticle
		print "is_lepton returns:" , tParticle.is_lepton()
		print "is_quark returns:" , tParticle.is_quark()
		print "is_meson returns:" , tParticle.is_meson()
		print "is_hadron returns:" , tParticle.is_hadron()
		print "is_baryon returns:" , tParticle.is_baryon()
		print "is_fermion returns:" , tParticle.is_fermion()
		print "is_gauge_boson returns:" , tParticle.is_gauge_boson()
		print "is_boson returns:" , tParticle.is_boson()
		print "is_anti_particle returns:" , tParticle.is_anti_particle()
		assertions.pause(__name__)
	print "\nFinished testing is__'' functions on all possible particles."
	assertions.pause(__name__)

	##Thoroughly test of which_particle_most_energic:##
	print "\n--------------------------------------------------\n"
	print "Testing thoroughly which_particle_most_energic:\n"
	print "Given two particles with: E1 = E2, E1 > E2 and then E2 > E1."
	tResults1, tResults2, tResults3 = [], [], []
	for ti in range(tNumTests):
		tRandE = random.random()
		tPT1 = particle(random.randint(1,6),fourVectors.fourVector(tRandE,0.0,0.0,0.0))
		tPT2 = particle(random.randint(1,6),fourVectors.fourVector(tRandE+1.0,0.0,0.0,0.0))
		tResults1.append(which_particle_most_energic(tPT1,tPT1))
		tResults2.append(which_particle_most_energic(tPT2,tPT1))
		tResults3.append(which_particle_most_energic(tPT1,tPT2))
	tFig = pyplot.figure()
	tAx = tFig.add_subplot(111)
	tTitle = r"$\rm{Testing\ which\_particle\_most\_energic\ for\ }$"
	tTitle += ("$"+str(tNumTests)+"$").encode('string-escape') + r"$\ \rm{iterations}$"
	pyplot.title(tTitle)
	pyplot.ylabel(r"$\rm{Average\ dipole\ index}$")
	pyplot.xlabel(r"$\rm{Iteration}$")
	pyplot.plot(tResults2, linewidth = 2, color = 'red', label = r"$\rm{E_{1} > E_{2}}$")
	pyplot.plot(tResults3, linewidth = 2, color = 'green', label = r"$\rm{E_{2} > E_{1}}$")
	pyplot.axhline(sum(tResults1)/float(tNumTests), linewidth = 2, color = 'blue')
	tBlueLine = lines.Line2D([], [], color='blue', linewidth = 2, label = r"$\rm{E_{1} = E_{2}}$")
	pyplot.ylim(-0.1,1.1)
	tTheHandles, tTheLabels = tAx.get_legend_handles_labels()
	tMyHandles = [tBlueLine] + tTheHandles
	tMyLabels = [r"$\rm{E_{1} = E_{2}}$"] + tTheLabels
	pyplot.legend(tMyHandles,tMyLabels)
	assertions.show_graph()
	print "The average for E1 = E2 is:", sum(tResults1)/float(tNumTests)
	print "This should be approximately 0.5?"
	assertions.pause(__name__)
	print "The average for E1 > E2 is:", sum(tResults2)/float(tNumTests)
	print "The average for E2 > E1 is:", sum(tResults3)/float(tNumTests)
	if ((sum(tResults2)/float(tNumTests) == 0) and (sum(tResults3)/float(tNumTests) == 1)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing thoroughly which_particle_most_energic."
	assertions.pause(__name__)

	##Thoroughly test of which_particle_to_recoil_g_emission:##
	print "\n--------------------------------------------------\n"
	print "Testing thoroughly which_particle_to_recoil_g_emission:\n"
	tResults = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
	for ti in range(tNumTests):
		tRandE = random.random()
		tG1 = particle(21,fourVectors.fourVector(tRandE,0.0,0.0,0.0))
		tG2 = particle(21,fourVectors.fourVector(tRandE+1.0,0.0,0.0,0.0))
		tResults[0].append(which_particle_to_recoil_g_emission(tG1,tG1))
		tResults[1].append(which_particle_to_recoil_g_emission(tG2,tG1))
		tResults[2].append(which_particle_to_recoil_g_emission(tG1,tG2))
		tG3 = particle(21,fourVectors.fourVector(tRandE,0.0,0.0,0.0))
		tQBar1 = particle(-random.randint(1,6),fourVectors.fourVector(tRandE*random.randrange(0.0,2.0),0.0,0.0,0.0)) #Energy randomly greater or less than q.
		tResults[3].append(which_particle_to_recoil_g_emission(tG3,tQBar1))
		tResults[4].append(which_particle_to_recoil_g_emission(tQBar1,tG3))
		tQ1 = particle(random.randint(1,6),fourVectors.fourVector(tRandE*random.randrange(0.0,2.0),0.0,0.0,0.0)) #Energy randomly greater or less than q.
		tResults[5].append(which_particle_to_recoil_g_emission(tG3,tQ1))
		tResults[6].append(which_particle_to_recoil_g_emission(tQ1,tG3))
		tQ1 = particle(random.randint(1,6),fourVectors.fourVector(tRandE,0.0,0.0,0.0))
		tQBar1 = particle(-random.randint(1,6),fourVectors.fourVector(tRandE,0.0,0.0,0.0))
		tQ3 = particle(random.randint(1,6),fourVectors.fourVector(tRandE*2.0,0.0,0.0,0.0))
		tQBar3 = particle(-random.randint(1,6),fourVectors.fourVector(tRandE*2.0,0.0,0.0,0.0))
		tResults[7].append(which_particle_to_recoil_g_emission(tQBar1,tQ1)) ##Same energy.
		tResults[8].append(which_particle_to_recoil_g_emission(tQBar3,tQ1))
		tResults[9].append(which_particle_to_recoil_g_emission(tQBar1,tQ3))
		tResults[10].append(which_particle_to_recoil_g_emission(tQ1,tQBar1)) ##Same energy.
		tResults[11].append(which_particle_to_recoil_g_emission(tQ3,tQBar1))
		tResults[12].append(which_particle_to_recoil_g_emission(tQ1,tQBar3))
	tFig = pyplot.figure()
	tAx = tFig.add_subplot(111)
	tTitle = r"$\rm{Testing\ which\_particle\_to\_recoil\ with\ g-g\ for\ }$"
	tTitle += ("$"+str(tNumTests)+"$").encode('string-escape') + r"$\ \rm{iterations}$"
	pyplot.title(tTitle)
	pyplot.ylabel(r"$\rm{Average\ dipole\ index}$")
	pyplot.xlabel(r"$\rm{Iteration}$")
	pyplot.plot(tResults[1], linewidth = 2, color = 'red', label = r"$\rm{E_{1} > E_{2}}$")
	pyplot.plot(tResults[2], linewidth = 2, color = 'green', label = r"$\rm{E_{2} > E_{1}}$")
	pyplot.axhline(sum(tResults[0])/float(tNumTests), linewidth = 2, color = 'blue')
	tBlueLine = lines.Line2D([], [], color='blue', linewidth = 2, label = r"$\rm{E_{1} = E_{2}}$")
	pyplot.ylim(-0.1,1.1)
	tTheHandles, theLabels = tAx.get_legend_handles_labels()
	tMyHandles = [tBlueLine] + tTheHandles
	tMyLabels = [r"$\rm{E_{1} = E_{2}}$"] + theLabels
	pyplot.legend(tMyHandles,tMyLabels)
	assertions.show_graph()
	print "For the gluons:"
	print "The average for E1 = E2 is:", sum(tResults[0])/float(tNumTests)
	print "This should be approximately 0.5?"
	assertions.pause(__name__)
	print "The average for E1 > E2 is:", sum(tResults[1])/float(tNumTests)
	print "The average for E2 > E1 is:", sum(tResults[2])/float(tNumTests)
	if ((sum(tResults[1])/float(tNumTests) == 0) and (sum(tResults[2])/float(tNumTests) == 1)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	assertions.pause(__name__)
	pyplot.figure()
	tTitle = r"$\rm{Testing\ which\_particle\_to\_recoil\ with\ g-\bar{q}\ for\ }$"
	tTitle += ("$"+str(tNumTests)+"$").encode('string-escape') + r"$\ \rm{iterations}$"
	pyplot.title(tTitle)
	pyplot.ylabel(r"$\rm{Average\ dipole\ index}$")
	pyplot.xlabel(r"$\rm{Iteration}$")
	pyplot.plot(tResults[3], linewidth = 2, color = 'red', label = r"$\rm{g-\bar{q}}$")
	pyplot.plot(tResults[4], linewidth = 2, color = 'green', label = r"$\rm{\bar{q}-g}$")
	pyplot.ylim(-0.1,1.1)
	pyplot.legend()
	assertions.show_graph()
	pyplot.figure()
	tTitle = r"$\rm{Testing\ which\_particle\_to\_recoil\ with\ g-q\ for\ }$"
	tTitle += ("$"+str(tNumTests)+"$").encode('string-escape') + r"$\ \rm{iterations}$"
	pyplot.title(tTitle)
	pyplot.ylabel(r"$\rm{Average\ dipole\ index}$")
	pyplot.xlabel(r"$\rm{Iteration}$")
	pyplot.plot(tResults[5], linewidth = 2, color = 'red', label = r"$\rm{g-q}$")
	pyplot.plot(tResults[6], linewidth = 2, color = 'green', label = r"$\rm{q-g}$")
	pyplot.ylim(-0.1,1.1)
	pyplot.legend()
	assertions.show_graph()
	print "For any g-q/qBar or q/qBar-g pair the gluon should recoil."
	print "For the gluon and anti-quark:"
	print "The average of the g - a-q test is:", sum(tResults[3])/float(tNumTests)
	print "The average of the a-q - g test is:", sum(tResults[4])/float(tNumTests)
	print "For the gluon and quark:"
	print "The average of the g - q test is:", sum(tResults[5])/float(tNumTests)
	print "The average of the q - g test is:", sum(tResults[6])/float(tNumTests)
	tCheck1 = ((sum(tResults[3])/float(tNumTests) == 0) and (sum(tResults[4])/float(tNumTests) == 1))
	tCheck2 = ((sum(tResults[5])/float(tNumTests) == 0) and (sum(tResults[6])/float(tNumTests) == 1))
	if (tCheck1 and tCheck2):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	assertions.pause(__name__)
	
	print "For qBar-q and q-qBar:"
	pyplot.figure()
	tTitle = r"$\rm{Testing\ which\_particle\_to\_recoil\ with\ \bar{q}-q\ for\ }$"
	tTitle += ("$"+str(tNumTests)+"$").encode('string-escape') + r"$\ \rm{iterations}$"
	pyplot.title(tTitle)
	pyplot.ylabel(r"$\rm{Average\ dipole\ index}$")
	pyplot.xlabel(r"$\rm{Iteration}$")
	pyplot.axhline(sum(tResults[7])/float(tNumTests), linewidth = 2, color = 'blue')
	pyplot.axhline(sum(tResults[8])/float(tNumTests), linewidth = 2, color = 'red')
	pyplot.axhline(sum(tResults[9])/float(tNumTests), linewidth = 2, color = 'green')
	tBlueLine = lines.Line2D([], [], color='blue', linewidth = 2, label = r"$\rm{E_{\bar{q}} = E_{q}}$")
	tRedLine = lines.Line2D([], [], color='red', linewidth = 2, label = r"$\rm{E_{\bar{q}} > E_{q}}$")
	tGreenLine = lines.Line2D([], [], color='green', linewidth = 2, label = r"$\rm{E_{q} > E_{\bar{q}}}$")
	pyplot.ylim(-0.1,1.1)
	pyplot.legend([tBlueLine,tRedLine,tGreenLine], [r"$\rm{E_{\bar{q}} = E_{q}}$",r"$\rm{E_{\bar{q}} > E_{q}}$",r"$\rm{E_{q} > E_{\bar{q}}}$"])
	assertions.show_graph()
	pyplot.figure()
	tTitle = r"$\rm{Testing\ which\_particle\_to\_recoil\ with\ q-\bar{q}\ for\ }$"
	tTitle += ("$"+str(tNumTests)+"$").encode('string-escape') + r"$\ \rm{iterations}$"
	pyplot.title(tTitle)
	pyplot.ylabel(r"$\rm{Average\ dipole\ index}$")
	pyplot.xlabel(r"$\rm{Iteration}$")
	pyplot.axhline(sum(tResults[10])/float(tNumTests), linewidth = 2, color = 'blue')
	pyplot.axhline(sum(tResults[11])/float(tNumTests), linewidth = 2, color = 'red')
	pyplot.axhline(sum(tResults[12])/float(tNumTests), linewidth = 2, color = 'green')
	tBlueLine = lines.Line2D([], [], color='blue', linewidth = 2, label = r"$\rm{E_{q} = E_{\bar{q}}}$")
	tRedLine = lines.Line2D([], [], color='red', linewidth = 2, label = r"$\rm{E_{q} > E_{\bar{q}}}$")
	tGreenLine = lines.Line2D([], [], color='green', linewidth = 2, label = r"$\rm{E_{\bar{q}} > E_{q}}$")
	pyplot.ylim(-0.1,1.1)
	pyplot.legend([tBlueLine,tRedLine,tGreenLine], [r"$\rm{E_{q} = E_{\bar{q}}}$",r"$\rm{E_{q} > E_{\bar{q}}}$",r"$\rm{E_{\bar{q}} > E_{q}}$"])
	assertions.show_graph()
	print "The average of E1 = E2 for qBar-q is:", sum(tResults[7])/float(tNumTests)
	print "The average of E1 = E2 for q-qBar is:", sum(tResults[10])/float(tNumTests)
	print "These should be approximately 0.5?"
	assertions.pause(__name__)
	print "The average of E1 > E2 for qBar-q is:", sum(tResults[8])/float(tNumTests)
	print "The average of E2 > E1 for qBar-q is:", sum(tResults[9])/float(tNumTests)
	print "These should be approximately 0.8 and 0.2?"
	print "The average of E1 > E2 for q-qBar is:", sum(tResults[11])/float(tNumTests)
	print "The average of E2 > E1 for q-qBar is:", sum(tResults[12])/float(tNumTests)
	print "These should be approximately 0.8 and 0.2?"
	print "\nFinished testing thoroughly which_particle_to_recoil_g_emission."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "///////////////////////////////////"
	print "Finished checking particles module!"
	print "///////////////////////////////////"