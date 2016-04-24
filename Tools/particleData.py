####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module containing classes for storing particle data and cataloguing particles by PDG values."""

##Import required modules:##
import assertions

print "\n////////////////////////////"
print "Loading particleData module:"
print "////////////////////////////\n"

##Functions:##

def check_is_particleDataEntry(toCheck):
	"""A function to check for an instance of the particleDataEntry class."""
	return isinstance(toCheck,particleDataEntry)

def check_is_allParticles(toCheck):
	"""A function to check for an instance of the allParticles class."""
	return isinstance(toCheck,allParticles)

def get_short_particle_name(pCodeOrLName):
	"""A function to return a short name for a particle given its code or long name."""
	assert ((type(pCodeOrLName) == str) or (type(pCodeOrLName) == int))
	if (type(pCodeOrLName) == str):
		__longName = pCodeOrLName
	else:
		__longName = knownParticles.get_name_from_code(pCodeOrLName)
	if (len(__longName) <= 4):
		__shortName = __longName
	else:
		if (__longName[0:5] == "anti-"):
			if ("muon" in __longName):
				__shortName = "a-mu" ##Particular case differs from most.
			else:
				__shortName = "a-" + __longName[5:8]
		elif ("neutrino" in __longName):
			__shortName = __longName[0] + "-nu"
		else:
			__shortName = __longName[0:3]
	return __shortName

##Classes:##

class particleDataEntry(object):
	"""A class to handle particle properties (name,mass,width,charge,spin) of a PDG entry.."""

	def __init__(self,nameToSet,massToSet = 0.0,widthToSet = 0.0,chargeToSet = 0.0,spinToSet = 0.0):
		"""A function to initiate a PDG entry. Only particles not anti-knownParticles."""
		assert (type(nameToSet) == str)
		assert assertions.all_are_numbers([massToSet,widthToSet,chargeToSet,spinToSet])
		self.__name  = nameToSet
		self.__mass  = assertions.force_float_number(massToSet)
		self.__width = assertions.force_float_number(widthToSet)
		self.__charge = assertions.force_float_number(chargeToSet)
		self.__spin = assertions.force_float_number(spinToSet)

	def get_name(self):
		"""A function to return the name of a PDG entry."""
		return self.__name

	def get_mass(self):
		"""A function to return the mass of a PDG entry."""
		return self.__mass

	def get_width(self):
		"""A function to return the width of a PDG entry."""
		return self.__width

	def get_charge(self):
		"""A function to return the charge of a PDG entry."""
		return self.__charge

	def get_spin(self):
		"""A function to return the spin of a PDG entry."""
		return self.__spin

class allParticles(object):
	"""A class containing PDG entries: a dictionary of particles with their names and properties."""
	##Only stores particles but the functions handle anti-particles called with -ve PDG code.
	##Natural units are used and masses are expressed in units of GeV: Width in GeV.
	##Values source: PDG online data sheets, 5/11/15.

	def __init__(self):
		"""A function to initiate a dictionary of particles with their names and properties."""
		self.__pDEs = {
			1:	particleDataEntry('d-quark',0.0048,0.0,-1.0/3.0,0.5),
			2:	particleDataEntry('u-quark',0.0023,0.0,2.0/3.0,0.5),
			3:	particleDataEntry('s-quark',0.095,0.0,-1.0/3.0,0.5),
			4:	particleDataEntry('c-quark',1.275,0.0,2.0/3.0,0.5),
			5:	particleDataEntry('b-quark',4.66,0.0,-1.0/3.0,0.5),		##Using the 1S decay mass.
			6:	particleDataEntry('t-quark',173.21,1.41,2.0/3.0,0.5),
			11: particleDataEntry('electron',0.000511,0.0,-1.0,0.5),
            12:  particleDataEntry('electron-neutrino',0.0,0.0,0.0,0.5),
            13:  particleDataEntry('muon',0.106,0.0,-1.0,0.5),
            14:  particleDataEntry('muon-neutrino',0.0,0.0,0.0,0.5),
            15:  particleDataEntry('tau',1.777,2.36e-12,-1.0,0.5),   ##Width uncertain.
            16:  particleDataEntry('tau-neutrino',0.0,0.0,0.0,0.5),
			21:	particleDataEntry('gluon',0.0,0.0,0.0,1.0),
			22:	particleDataEntry('photon',0.0,0.0,0.0,1.0),
			23: particleDataEntry('Z-boson',91.188,2.49,0.0,1.0)}

		self.__quarks   = [1,2,3,4,5,6]
		self.__leptons  = [11,12,13,14,15,16]
		self.__mesons   = []
		self.__baryons  = []
		self.__gaugeBosons = [21,22,23]
		self.__hadrons  = self.__mesons + self.__baryons
		self.__bosons = self.__mesons + self.__gaugeBosons
		self.__fermions = self.__quarks + self.__leptons + self.__baryons
		self.__allCodes = self.__fermions + self.__bosons
		self.__hadrons.sort(), self.__bosons.sort() , self.__fermions.sort(), self.__allCodes.sort()

	def has_key(self,lookUpCode):
		"""A function to check if the particles class contains a given PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		if self.__pDEs.has_key(abs(lookUpCode)):
			return True
		return False

	def has_name(self,lookUpName):
		"""A function to check if the particles class contains a given name."""
		assert (type(lookUpName) == str)
		for __lookUpCode in self.__pDEs:
			if (self.__pDEs[__lookUpCode].get_name() == lookUpName):
				return True
		return False

	def get_name_from_code(self,lookUpCode):
		"""A function to return the name of a particle given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		if (lookUpCode > 0.0):
			return self.__pDEs[lookUpCode].get_name()
		return "anti-" + self.__pDEs[abs(lookUpCode)].get_name()

	def get_mass_from_code(self,lookUpCode):
		"""A function to return the mass of a particle given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return self.__pDEs[abs(lookUpCode)].get_mass()


	def get_width_from_code(self,lookUpCode):
		"""A function to return the width of a particle given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return self.__pDEs[abs(lookUpCode)].get_width()

	def get_charge_from_code(self,lookUpCode):
		"""A function to return the charge of a particle given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		if (lookUpCode > 0.0):
			return self.__pDEs[lookUpCode].get_charge()
		return (-1.0 * self.__pDEs[abs(lookUpCode)].get_charge())

	def get_spin_from_code(self,lookUpCode):
		"""A function to return the spin of a particle given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return self.__pDEs[abs(lookUpCode)].get_spin()

	def get_code_from_name(self,lookUpName):
		"""A function to return the PDG code of a particle given its name."""
		assert (type(lookUpName) == str)
		__nameExists = False
		if (lookUpName[0:5] == "anti-"):
			for __PDGCode in self.__pDEs:
				if (self.__pDEs[__PDGCode].get_name() == lookUpName[5:]):
					__nameExists = True
					return -1 * __PDGCode
		else:
			for __PDGCode in self.__pDEs:
				if (self.__pDEs[__PDGCode].get_name() == lookUpName):
					__nameExists = True
					return __PDGCode
		assert __nameExists

	def is_lepton(self,lookUpCode):
		"""A function to check if a particle is a lepton given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
	 	return (self.__leptons.count(abs(lookUpCode)) > 0.0)

	def is_quark(self,lookUpCode):
		"""A function to check if a particle is a quark given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__quarks.count(abs(lookUpCode)) > 0.0)

	def is_meson(self,lookUpCode):
		"""A function to check if a particle is a meson given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__mesons.count(abs(lookUpCode)) > 0.0)

	def is_hadron(self,lookUpCode):
		"""A function to check if a particle is a hadron given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__hadrons.count(abs(lookUpCode)) > 0.0)

	def is_baryon(self,lookUpCode):
		"""A function to check if a particle is a baryon given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__baryons.count(abs(lookUpCode)) > 0.0)

	def is_fermion(self,lookUpCode):
		"""A function to check if a particle is a fermion given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__fermions.count(abs(lookUpCode)) > 0.0)

	def is_gauge_boson(self,lookUpCode):
		"""A function to check if a particle is a gauge boson given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__gaugeBosons.count(abs(lookUpCode)) > 0.0)

	def is_boson(self,lookUpCode):
		"""A function to check if a particle is a boson given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		assert self.has_key(abs(lookUpCode))
		return (self.__bosons.count(abs(lookUpCode)) > 0.0)

	@staticmethod
	def is_anti_particle(lookUpCode):
		"""A function to check if a particle is an anti-particle given its PDG code."""
		assert assertions.all_are_numbers([lookUpCode])
		return (lookUpCode < 0.0)

	def get_known_quarks(self):
		"""A function to return a list of known quark PDG codes."""
		return self.__quarks

	def get_known_fermions(self):
		"""A function to return a list of known fermion PDG codes."""
		return self.__fermions

	def get_all_codes(self):
		"""A function to return a list of all known PDG codes."""
		return self.__allCodes

##Create an instance of allParticles() to be used whenever the module is included:
knownParticles = allParticles()

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////////"
	print "Testing particleData module:"
	print "////////////////////////////"
	assertions.pause(__name__)

	##Setup here:##
	print "\nGenerating test values..."
	tPDE1 = particleDataEntry('testParticle',10.5,2.5,1.0/3.0,1.0/2.0)
	tString = "only a string"
	tShortNameList = [1,-1,2,-2,3,-3,4,-4,5,-5,6,-6,11,-11,12,13,-13,14,15,-15,16,21,22,23]

	##Test check_is_'' functions:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_'' functions:\n"
	print "Calling check_is_particleDataEntry on particleDataEntry returns:", check_is_particleDataEntry(tPDE1)
	print "Calling check_is_particleDataEntry on a string returns:", check_is_particleDataEntry(tString)
	print "Calling check_is_particleDataEntry on a float returns:", check_is_particleDataEntry(2.005)
	print "Calling check_is_allParticles on allParticles returns:", check_is_allParticles(knownParticles)
	print "Calling check_is_allParticles on a string returns:", check_is_allParticles(tString)
	print "Calling check_is_allParticles on a float returns:", check_is_allParticles(2.005)
	tResults1 = [check_is_particleDataEntry(tPDE1),check_is_particleDataEntry(tString),check_is_particleDataEntry(2.005)]
	tResults2 = [check_is_allParticles(knownParticles),check_is_allParticles(tString),check_is_allParticles(2.005)]
	tResults = tResults1 + tResults2
	if ((sum(tResults) == 2) and (tResults[0] == True) and (tResults[3] == True)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_'' functions."
	assertions.pause(__name__)

	##Test get_short_particle_name function:##
	print "\n--------------------------------------------------\n"
	print "Testing get_short_particle_name function:\n"
	for tCode in tShortNameList:
		print "Calling on:", tCode, "\ngives:", get_short_particle_name(tCode)
		testLongName = knownParticles.get_name_from_code(tCode)
		print "Calling via it's full name:", testLongName, "\ngives:", get_short_particle_name(testLongName)
		assertions.pause(__name__)
	print "\nFinished testing get_short_particle_name function."
	assertions.pause(__name__)

	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////////////"
	print "Testing particleDataEntry class:"
	print "////////////////////////////////"
	assertions.pause(__name__)

	##Test __init__ and get_'' functions of particleDataEntry class:##
	print "\n--------------------------------------------------\n"
	print "Testing __init__ and get_'' functions of particleDataEntry class:\n"
	print "Entry 1 has values:" , tPDE1.get_name(), tPDE1.get_mass(), tPDE1.get_width()
	print  "and:", tPDE1.get_charge(), tPDE1.get_spin()
	tResults = [tPDE1.get_name(), tPDE1.get_mass(), tPDE1.get_width(),tPDE1.get_charge(), tPDE1.get_spin()]
	if (tResults == ['testParticle',10.5,2.5,1.0/3.0,1.0/2.0]):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing __init__ and get_'' functions of particleDataEntry class."
	assertions.pause(__name__)

	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////////"
	print "Testing allParticles class:"
	print "///////////////////////////"
	assertions.pause(__name__)

	##__init__ tested implicitly.

	##Test has_'' functions of allParticles class:##
	print "\n--------------------------------------------------\n"
	print "Testing has_'' functions of allParticles class:\n"
	print "Given the key 4 which exists returns :" , knownParticles.has_key(4)
	print "Given the key -4 returns :" , knownParticles.has_key(-4)
	print "Given the key 100 which doesn't exist returns :" , knownParticles.has_key(100)
	print "Given the key -100 returns :" , knownParticles.has_key(-100)
	print "Given the name 'photon' which exists returns :" , knownParticles.has_name('photon')
	print "Given the name 's-quark' which exists returns :" , knownParticles.has_name('s-quark')
	print "Given the name 'squark' which doesn't exist returns :" , knownParticles.has_name('squark')
	print "Given the name 'random' which doesn't exist returns :" , knownParticles.has_name('random')
	print "Given the name '' which doesn't exist returns :" , knownParticles.has_name('')
	tResults = [knownParticles.has_key(4),knownParticles.has_key(-4),knownParticles.has_key(100),knownParticles.has_key(-100)]
	tResults += [knownParticles.has_name('photon'),knownParticles.has_name('s-quark'),knownParticles.has_name('squark')]
	tResults += [knownParticles.has_name('random'),knownParticles.has_name('')]
	if (sum(tResults) == 4):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing has_'' functions of allParticles class."
	assertions.pause(__name__)

	##Test get_'' functions of allParticles class:##
	print "\n--------------------------------------------------\n"
	print "Testing get_'' functions of allParticles class:\n"
	print "Using the s-quark which has code 3:"
	print "Calling get_name_from_code returns:" , knownParticles.get_name_from_code(3)
	print "Calling get_name_from_code with -ve returns:" , knownParticles.get_name_from_code(-3)
	print "Calling get_mass_from_code returns:" , knownParticles.get_mass_from_code(3)
	print "Calling get_mass_from_code with -ve returns:" , knownParticles.get_mass_from_code(-3)
	print "Calling get_width_from_code returns:" , knownParticles.get_width_from_code(3)
	print "Calling get_width_from_code with -ve returns:" , knownParticles.get_width_from_code(-3)
	tResults = [knownParticles.get_name_from_code(3), knownParticles.get_name_from_code(-3), knownParticles.get_mass_from_code(3)]
	tResults += [knownParticles.get_mass_from_code(-3),knownParticles.get_width_from_code(3),knownParticles.get_width_from_code(-3)]
	assertions.pause(__name__)
	print "\nCalling get_charge_from_code returns:" , knownParticles.get_charge_from_code(3)
	print "Calling get_charge_from_code with -ve returns:" , knownParticles.get_charge_from_code(-3)
	print "Calling get_spin_from_code returns:" , knownParticles.get_spin_from_code(3)
	print "Calling get_spin_from_code with -ve returns:" , knownParticles.get_spin_from_code(-3)
	print "Calling get_code_from_name returns:" , knownParticles.get_code_from_name("s-quark")
	print "Calling get_code_from_name with 'anti-' returns:" , knownParticles.get_code_from_name("anti-s-quark")
	tResults += [knownParticles.get_charge_from_code(3),knownParticles.get_charge_from_code(-3),knownParticles.get_spin_from_code(3)]
	tResults += [knownParticles.get_spin_from_code(-3),knownParticles.get_code_from_name("s-quark"),knownParticles.get_code_from_name("anti-s-quark")]
	if (tResults == ['s-quark','anti-s-quark',0.095,0.095,0.0,0.0,-1.0/3.0,1.0/3.0,0.5,0.5,3,-3]):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_'' functions of allParticles class."
	assertions.pause(__name__)

	##Test is_'' functions of allParticles class:##
	print "\n--------------------------------------------------\n"
	print "Testing is_'' functions of allParticles class:\n"
	print "Using the s-quark which has code 3:"
	print "Calling is_lepton returns:" , knownParticles.is_lepton(3)
	print "Calling is_lepton with -ve returns:" , knownParticles.is_lepton(-3)
	print "Calling is_quark returns:" , knownParticles.is_quark(3)
	print "Calling is_quark with -ve returns:" , knownParticles.is_quark(-3)
	print "Calling is_meson returns:" , knownParticles.is_meson(3)
	print "Calling is_meson with -ve returns:" , knownParticles.is_meson(-3)
	tResults = [knownParticles.is_lepton(3),knownParticles.is_lepton(-3),knownParticles.is_quark(3)]
	tResults += [knownParticles.is_quark(-3),knownParticles.is_meson(3),knownParticles.is_meson(-3)]
	assertions.pause(__name__)
	print "\nCalling is_hadron returns:" , knownParticles.is_hadron(3)
	print "Calling is_hadron with -ve returns:" , knownParticles.is_hadron(-3)
	print "Calling is_baryon returns:" , knownParticles.is_baryon(3)
	print "Calling is_baryon with -ve returns:" , knownParticles.is_baryon(-3)
	print "Calling is_fermion returns:" , knownParticles.is_fermion(3)
	print "Calling is_fermion with -ve returns:" , knownParticles.is_fermion(-3)
	tResults += [knownParticles.is_hadron(3),knownParticles.is_hadron(-3),knownParticles.is_baryon(3)]
	tResults += [knownParticles.is_baryon(-3),knownParticles.is_fermion(3),knownParticles.is_fermion(-3)]
	assertions.pause(__name__)
	print "\nCalling is_gauge_boson returns:" , knownParticles.is_gauge_boson(3)
	print "Calling is_gauge_boson with -ve returns:" , knownParticles.is_gauge_boson(-3)
	print "Calling is_boson returns:" , knownParticles.is_boson(3)
	print "Calling is_boson with -ve returns:" , knownParticles.is_boson(-3)
	print "Calling is_anti_particle returns:" , knownParticles.is_anti_particle(3)
	print "Calling is_anti_particle with -ve returns:" , knownParticles.is_anti_particle(-3)
	tResults += [knownParticles.is_gauge_boson(3),knownParticles.is_gauge_boson(-3),knownParticles.is_boson(3)]
	tResults += [knownParticles.is_boson(-3),knownParticles.is_anti_particle(3),knownParticles.is_anti_particle(-3)]
	if (sum(tResults) == 5):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing is_'' functions of particleData class."
	assertions.pause(__name__)

	##Test get_known_quarks, get_known_fermions and get_all_codes function of allParticles class:##
	print "\n--------------------------------------------------\n"
	print "Testing get_known_quarks, get_known_fermions and get_all_codes function of allParticles class:\n"
	print "Calling get_known_quarks returns:", knownParticles.get_known_quarks()
	print "Calling get_known_fermions returns:", knownParticles.get_known_fermions()
	print "Calling get_all_codes returns:", knownParticles.get_all_codes()
	tResult = ((len(knownParticles.get_known_quarks()) == 6) and (len(knownParticles.get_known_fermions()) == 12))
	if (tResult and (len(knownParticles.get_all_codes()) == 15)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing get_known_quarks, get_known_fermions and get_all_codes function of allParticles class."
	assertions.pause(__name__)
	
	##Done testing:##
	print "\n---------------------------------------------\n"
	print "//////////////////////////////////////"
	print "Finished checking particleData module!"
	print "//////////////////////////////////////"