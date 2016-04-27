####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for handling chains in dipole showering."""

##Import required modules:##
import random
import assertions
import constants
import counters
import dataLoggers
import particleData
import fourVectors
import particles
import kinematics
import sudakovs
import dipoles
import quarkPairs

print "\n//////////////////////"
print "Loading chains module:"
print "//////////////////////\n"

##Set up data loggers:
producedQuarkCodes = dataLoggers.dataLogger()

##Functions:##

def check_is_chain(toCheck):
	"""A function to check for an instance of the chain class."""
	return isinstance(toCheck,chain)

##Classes:##

class chain(object):
	"""A class to handle a chain."""
	##Dipoles are accessed by list index. The same particle can be in two neighbouring dipoles.
	##e.g. [[a,b],[b,c],[c,d],[d,e]], where for a closed chain e == a.

	def __init__(self,orderedListOfParticles,isLoop = False,maxPperpSquaredFromSplit = None):
		"""A function to initialise a chain using the ordered list of particles given."""
		##Accepting an ordered list of particles allows the initialisation of a chain of any length at any stage.
		for particle in orderedListOfParticles:
			assert particles.check_is_particle(particle)
		assert (isLoop or not(isLoop))
		assert ((maxPperpSquaredFromSplit == None) or (assertions.all_are_numbers([maxPperpSquaredFromSplit])))
		self.__chainList = []
		for i in range(len(orderedListOfParticles) - 1): ##Iterate through all particle pairs.
			self.__chainList.append(dipoles.dipole(orderedListOfParticles[i],orderedListOfParticles[i+1]))
		self.__isLoop = isLoop
		if (maxPperpSquaredFromSplit == None):
			if (len(self.__chainList) == 1):
				self.__maxPperpSquared = self.__chainList[0].get_mass_squared()
			else:
				print "\n! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
				print "Not sure what to do here yet!"
				#assert False
				self.__maxPperpSquared = self.__chainList[0].get_mass_squared()
				assertions.pause_loading_module()
				##Or do you find that of the entire chain. Or do you find that of each -> Need list?
		else:
			self.__maxPperpSquared = maxPperpSquaredFromSplit
		self.__cutOff = constants.cut_off_energy()*constants.cut_off_energy()
		self.__nextPperpSquared = None ##Initiation value.
		self.__showeringCompleted = False

	def set_showering_completed(self):
		"""A function to set the showering completed status to True."""
		self.__showeringCompleted = True

	def showering_completed(self):
		"""A function to check if the showering for the entire chain is completed."""
		return self.__showeringCompleted

	def check_is_closed(self):
		"""A function to check whether the chain is a closed loop or not."""
		return self.__isLoop

	def get_chain_index(self,indexIn):
		"""A function to check a chain index, taking into account if the chain is a loop."""
		assert (type(indexIn) == int)
		__chainLength = len(self.__chainList)
		__chainIndex = indexIn%__chainLength
		##__chainLength runs from 1 not 0 hence < not <= below.
		assert ((__chainIndex < __chainLength) and (__chainIndex >= 0))
		return __chainIndex

	def __getitem__(self,indexIn):
		"""A function to return one of the dipoles in the chain using its chain-list index."""
		 ##Don't need a __setitem__ as only set when making chain. Other changes are made to the list.
		assert (type(indexIn) == int)
		__listIndexOfDipole = self.get_chain_index(indexIn)
		return self.__chainList[__listIndexOfDipole]

	def get_chain_list(self):
		"""A function to return the chain list of dipoles."""
		return self.__chainList

	def get_next_Pperp_squared(self):
		"""A function to return the next value of PperpSquared currently being used by a chain."""
		return self.__nextPperpSquared

	def get_max_Pperp_squared(self):
		"""A function to return the maximum value of PperpSquared currently being used by a chain."""
		return self.__maxPperpSquared

	def __str__(self):
		"""A function to return a string of a chain object."""
		__stringOfChain = "\n[chain-"
		for dipole in self.__chainList:
			__stringOfChain += (str(dipole) + " |\n")
		__stringOfChain += ("   Loop: " + str(self.__isLoop) + ", Cut-off: " + str(self.__cutOff) + "\n-]")
		return __stringOfChain

	def simple_str(self):
		"""A function to return a simplified string of a chain object."""
		##Uses arrows at beginning and end to indicate a closed/looped chain.
		__simpleStringOfChain = "["
		for dipole in self.__chainList:
			__simpleStringOfChain += dipole.simple_str()
		__simpleStringOfChain += "]"
		if self.__isLoop:
			return "<--" + __simpleStringOfChain + "-->"
		else:
			return __simpleStringOfChain

	def __repr__(self):
		"""A function to return a string of a chain object."""
		__reprOfChain = "\n<chain(["
		for dipole in self.__chainList:
			__reprOfChain += (repr(dipole) + " |\n")
		__reprOfChain += ("   Loop: " + str(self.__isLoop) + ", Cut-off: " + str(self.__cutOff) + "])>")
		return __reprOfChain

	def run_sudakovs(self):
		"""A function to run the sudakov module for each dipole in the chain."""
		##If self.__maxPperpSquared > self.__cutOff, still run as process_sudakovs will catch.
		##sudakovList contains lists of Pperp^2, y and processCode for each dipole.
		self.__sudakovLists = [[],[],[]]
		for dipole in self.__chainList:
			__S123 = dipole.get_mass_squared()
			__code1, __code2 = dipole[0].get_code(), dipole[1].get_code()
			__PperpSquared,__y,__processCode = sudakovs.solve(self.__maxPperpSquared,__S123,__code1,__code2)
			self.__sudakovLists[0].append(__PperpSquared)
			self.__sudakovLists[1].append(__y)
			self.__sudakovLists[2].append(__processCode)
	
	def process_sudakovs(self):
		"""A function to process the list of sudakov results."""
		if sum(self.__sudakovLists[2]) == 0: ##Implies all returned values below the cutoff.
			self.__maxPperpSquared = self.__cutOff ##Set to assure no further events occur for this string.
			return False
		else:
			self.__nextPperpSquared = max(self.__sudakovLists[0]) ##Find maximum PperpSquared value.
			##Check for any other occurances.
			__indices = [i for i, x in enumerate(self.__sudakovLists[0]) if x == self.__nextPperpSquared]
			if (len(__indices) > 1):
				self.__nextIndex = __indices[0]
			else: ##Randomly choose if more than one had the same maxPperp^2.
				__randIndex = int(round(random.random()*len(__indices)-1)) ##-1 for index from length.
				self.__nextIndex = __indices[__randIndex] 
			self.__nextY = self.__sudakovLists[1][self.__nextIndex]
			self.__nextProcessCode = self.__sudakovLists[2][self.__nextIndex]
			self.__maxPperpSquared = self.__nextPperpSquared ##Set the new maximum to that currently occuring.
			return True
	
	def prepare_dipoles(self):
		"""A function to prepare the two relevant dipoles for updating."""
		assert (self.__maxPperpSquared > self.__cutOff)
		self.__chainList[self.__nextIndex].set_not_updated()
		if (self.__isLoop or (self.__nextIndex > 0)):
			self.__chainList[self.get_chain_index(self.__nextIndex - 1)].set_not_updated()
		elif (self.__isLoop or (self.__nextIndex < (len(self.__chainList) - 1))): ##-1 for index from length.
			self.__chainList[self.get_chain_index(self.__nextIndex + 1)].set_not_updated()
		self.__dipoleRecoilIndex = self.__chainList[self.__nextIndex].get_index_to_recoil(self.__nextProcessCode)
		if (self.__dipoleRecoilIndex == 0): ## LHS recoils.
			if self.__isLoop:
				self.__toRecoil = [self.get_chain_index(self.__nextIndex - 1),1] ##Store [dipole index, particle index].
				self.__chainList[self.get_chain_index(self.__nextIndex - 1)].set_not_updated()
			elif (self.__nextIndex == 0):
				self.__toRecoil = [None,None] ##Store [dipole index, particle index].
			else:
				self.__toRecoil = [self.__nextIndex - 1,1] ##Store [dipole index, particle index].
				self.__chainList[self.__nextIndex - 1].set_not_updated()
		else: ##RHS recoils.
			if self.__isLoop:
				self.__toRecoil = [self.get_chain_index(self.__nextIndex + 1),0]
				self.__chainList[self.get_chain_index(self.__nextIndex + 1)].set_not_updated()
			elif self.__nextIndex == (len(self.__chainList) - 1): ##-1 to get index from length.
				self.__toRecoil = [None,None]
			else:
				self.__toRecoil = [self.__nextIndex + 1,0]
				self.__chainList[self.__nextIndex + 1].set_not_updated()
	
	def run_kinematics(self):
		"""A function to run the kinematics module for a given chain event."""
		assert (self.__maxPperpSquared > self.__cutOff)
		self.__dipoleProduceIndex = (self.__dipoleRecoilIndex + 1)%2
		__v1 = self.__chainList[self.__nextIndex][self.__dipoleProduceIndex].get_four_momentum()
		__v3 = self.__chainList[self.__nextIndex][self.__dipoleRecoilIndex].get_four_momentum()
		__nV1, __nV2, __nV3 = kinematics.produce_new_vectors(__v1,__v3,self.__nextPperpSquared,self.__nextY)
		return __nV1, __nV2, __nV3

	def update_dipole_values(self,nPLHS,nDPLHS,nDPRHS,nPRHS):
		"""A function to update the values in the surrounding dipoles and insert the new one."""
		assert particles.check_is_particle(nPLHS)
		assert dipoles.check_is_dipole(nDPLHS)
		assert dipoles.check_is_dipole(nDPRHS)
		assert particles.check_is_particle(nPRHS)
		if ((not self.__isLoop) and (self.__nextIndex == 0)): ##No particle of dipole to left to update.
			pass
		else:
			__uLHIndex = self.get_chain_index(self.__nextIndex - 1)
			self.__chainList[__uLHIndex][1] = nPLHS.copy()
		if ((not self.__isLoop) and (self.__nextIndex == (len(self.__chainList) - 1))): ##-1 for list index from length.
			pass ##No particle of dipole to right to update.
		else:
			__uRHIndex = self.get_chain_index(self.__nextIndex + 1)
			self.__chainList[__uRHIndex][0] = nPRHS.copy()
		self.__chainList[self.__nextIndex] = nDPRHS.copy()
		self.__chainList.insert(self.__nextIndex,nDPLHS.copy())

	##~~Gluon Emission~~##

	def set_new_colours_g_prod(self,cC,p1,p2,p3):
		"""A function to set the new colours of the produced particles."""
		assert counters.check_is_counter(cC)
		for p in [p1,p2,p3]:
			assert particles.check_is_particle(p)
			assert p.__nonzero__()
		__kQs = particleData.knownParticles.get_known_quarks()
		__gC = particleData.knownParticles.get_code_from_name('gluon')
		if ((p1.get_code() in __kQs) and ((-p3.get_code()) in __kQs)): ##qqBar.
			__nC = cC.next()
			p2.set_colour(0,__nC)
			p2.set_colour(1,p1.get_colour(0))
			p3.set_colour(1,__nC)
		elif (((-p1.get_code()) in __kQs) and (p3.get_code() in __kQs)): ##qBarq.
			__nC = cC.next()
			p1.set_colour(1,__nC)
			p2.set_colour(0,__nC)
			p2.set_colour(1,p3.get_colour(0))
		elif ((p1.get_code() in __kQs) and (p3.get_code() == __gC)): ##qg.
			__nC = cC.next()
			p2.set_colour(0,__nC)
			p2.set_colour(1,p1.get_colour(0))
			p3.set_colour(1,__nC)
		elif (((-p1.get_code()) in __kQs) and (p3.get_code() == __gC)): ##qBarg.
			__nC = cC.next()
			p2.set_colour(0,p1.get_colour(1))
			p2.set_colour(1,__nC)
			p3.set_colour(0,__nC)
		elif ((p1.get_code() == __gC) and (p3.get_code() in __kQs)): ##gq.
			__nC = cC.next()
			p1.set_colour(1,__nC)
			p2.set_colour(0,__nC)
			p2.set_colour(1,p3.get_colour(0))
		elif ((p1.get_code() == __gC) and ((-p3.get_code()) in __kQs)): ##gqBar.
			__nC = cC.next()
			p1.set_colour(0,__nC)
			p2.set_colour(0,p3.get_colour(1))
			p2.set_colour(1,__nC)
		elif ((p1.get_code() == __gC) and (p3.get_code() == __gC)): ##gg.
			##Two possibilities so randomly choose which keeps colour code.
			__R = random.random()
			if (self.__dipoleRecoilIndex == 0): ##qBar at p3 end, i.e LHS.
				if (__R < 0.5): ##Choose option A.
					__nC = cC.next()
					p1.set_colour(1,__nC)
					p2.set_colour(0,__nC)
					p2.set_colour(1,p3.get_colour(0))
				else: ##Choose option B.
					__nC = cC.next()
					p2.set_colour(0,p1.get_colour(1))
					p2.set_colour(1,__nC)
					p3.set_colour(0,__nC)
			else: ##q at p3 end, i.e RHS.
				if (__R < 0.5): ##Choose option A.
					__nC = cC.next()
					p2.set_colour(0,__nC)
					p2.set_colour(1,p1.get_colour(0))
					p3.set_colour(1,__nC)
				else: ##Choose option B.
					__nC = cC.next()
					p1.set_colour(0,__nC)
					p2.set_colour(0,p3.get_colour(1))
					p2.set_colour(1,__nC)
		return p1, p2, p3

	def update_dipoles_g_prod_LHS_recoil(self,pC,cC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a gluon is produced with the LHS recoiling."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		__p1 = self.__chainList[self.__nextIndex][1].copy()
		__p2 = particles.particle('gluon',nV2,[0,0],[0,0],[0,0],1)
		__p3 = self.__chainList[self.__nextIndex][0].copy()
		__p1.set_four_momentum(nV1)
		__p3.set_four_momentum(nV3)
		__p2.set_unique_ID(pC.next()) ##Set unique particle ID when produced.
		__p2.set_produced_at(self.__maxPperpSquared) ##As now set to that of this produced particle.
		__p2.set_mother(0,__p1.get_unique_ID())
		__p2.set_mother(1,__p3.get_unique_ID())
		if __p1.get_child(0) != 0: ##i.e don't overwrite
			__1.set_child(0,__p2.get_unique_ID())
		if __p3.get_child(1) != 0: ##i.e don't overwrite
			__p3.set_child(1,__p2.get_unique_ID())
		__p1, __p2, __p3 = self.set_new_colours_g_prod(cC,__p1,__p2,__p3)
		__nDPLHS = self.__chainList[self.__nextIndex].copy()
		__nDPRHS = self.__chainList[self.__nextIndex].copy()
		__nDPLHS[0], __nDPLHS[1] = __p3.copy(), __p2.copy()
		__nDPRHS[0], __nDPRHS[1] = __p2.copy(), __p1.copy()
		self.update_dipole_values(__p3.copy(),__nDPLHS,__nDPRHS,__p1.copy())

	def update_dipoles_g_prod_RHS_recoil(self,pC,cC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a gluon is produced with the RHS recoiling."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		__p1 = self.__chainList[self.__nextIndex][0].copy()
		__p2 = particles.particle('gluon',nV2,[0,0],[0,0],[0,0],1)
		__p3 = self.__chainList[self.__nextIndex][1].copy()
		__p1.set_four_momentum(nV1)
		__p3.set_four_momentum(nV3)
		__p2.set_unique_ID(pC.next()) ##Set unique particle ID when produced.
		__p2.set_produced_at(self.__maxPperpSquared) ##As now set to that of this produced particle.
		__p2.set_mother(0,__p1.get_unique_ID())
		__p2.set_mother(1,__p3.get_unique_ID())
		if __p1.get_child(1) != 0: ##i.e don't overwrite
			__p1.set_child(1,__p2.get_unique_ID())
		if __p3.get_child(0) != 0: ##i.e don't overwrite
			__p3.set_child(0,__p2.get_unique_ID())
		__p1, __p2, __p3 = self.set_new_colours_g_prod(cC,__p1,__p2,__p3)
		__nDPLHS = self.__chainList[self.__nextIndex].copy()
		__nDPRHS = self.__chainList[self.__nextIndex].copy()
		__nDPLHS[0], __nDPLHS[1] = __p1.copy(), __p2.copy()
		__nDPRHS[0], __nDPRHS[1] = __p2.copy(), __p3.copy()
		self.update_dipole_values(__p1.copy(),__nDPLHS,__nDPRHS,__p3.copy())

	def update_dipoles_g_prod(self,pC,cC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a gluon is produced."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		if (self.__dipoleRecoilIndex == 0): ##LHS recoiling.
			self.update_dipoles_g_prod_LHS_recoil(pC,cC,nV1,nV2,nV3)
		else: ##RHS recoiling.
			self.update_dipoles_g_prod_RHS_recoil(pC,cC,nV1,nV2,nV3)

	##~~Gluon Split~~##

	def set_new_colours_g_split(self,cC,p1,p2,p3,p1C0B,p1C1B,hS):
		"""A function to set the new colours of the produced particles."""
		assert counters.check_is_counter(cC)
		for p in [p1,p2,p3]:
			assert particles.check_is_particle(p)
			assert p.__nonzero__()
		assert ((type(p1C0B) == int) and (type(p1C1B) == int))
		assert ((hS == "LHS") or (hS == "RHS"))
		__kQs = particleData.knownParticles.get_known_quarks()
		__gC = particleData.knownParticles.get_code_from_name('gluon')
		if (hS == "LHS"):
			if (p1.get_code() > 0): ##Q
				p1.set_colour(0,p1C0B)
				p2.set_colour(1,p1C1B)
			else: ##QBar
				p1.set_colour(0,p1C0B)
				p2.set_colour(1,p1C1B)
		elif (hS == "RHS"):
			if (p1.get_code() > 0): ##Q
				p1.set_colour(1,p1C0B)
				p2.set_colour(0,p1C1B)
			else: ##QBar
				p1.set_colour(1,p1C1B)
				p2.set_colour(0,p1C0B)
		return p1, p2, p3

	def update_dipoles_g_split_LHS_recoil(self,pC,cC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a gluon is split with the LHS recoiling."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		__p1B = self.__chainList[self.__nextIndex][1].copy()
		__p1BMs, __p1BCs = [__p1B.get_mother(0),__p1B.get_mother(1)], [__p1B.get_child(0),__p1B.get_child(1)]
		__qCode = quarkPairs.get_quark_code(self.__maxPperpSquared,self.__activeQCodes) ##As maxPperpSquared now set to S123 of splitting.
		producedQuarkCodes.store(__qCode)
		__p1C0B, __p1C1B = __p1B.get_colour(0), __p1B.get_colour(1)
		__p1 = particles.particle(__qCode,__p1B.get_four_momentum(),__p1BMs,__p1BCs,[0,0],__p1B.get_status_code())
		__p2 = particles.particle(-__qCode,nV2,[0,0],[0,0],[0,0],1)
		__p3 = self.__chainList[self.__nextIndex][0].copy()
		__p1.set_four_momentum(nV1)
		__p3.set_four_momentum(nV3)
		__p1.set_unique_ID(pC.next()) ##Set unique particle IDs when produced.
		__p2.set_unique_ID(pC.next())
		__p1.set_produced_at(self.__maxPperpSquared) ##As now set to that of this produced particle.
		__p2.set_produced_at(self.__maxPperpSquared)
		__p1.set_mother(0,__p3.get_unique_ID())
		__p2.set_mother(1,__p3.get_unique_ID())
		if __p3.get_child(1) != 0: ##i.e don't overwrite
			__p3.set_child(1,__p2.get_unique_ID())
		__p1, __p2, __p3 = self.set_new_colours_g_split(cC,__p1,__p2,__p3, __p1C0B, __p1C1B,"LHS")
		__nDPLHS = self.__chainList[self.__nextIndex].copy()
		__nDPRHS = self.__chainList[self.__nextIndex].copy()
		__nDPLHS[0], __nDPLHS[1] = __p3.copy(), __p2.copy()
		__nDPRHS[0], __nDPRHS[1] = __p2.copy(), __p1.copy()
		self.update_dipole_values(__p3.copy(),__nDPLHS,__nDPRHS,__p1.copy())
		self.__sideOfSplit = "LHS"

	def update_dipoles_g_split_RHS_recoil(self,pC,cC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a gluon is split with the RHS recoiling."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		__p1B = self.__chainList[self.__nextIndex][0].copy()
		__p1BMs, __p1BCs = [__p1B.get_mother(0),__p1B.get_mother(1)], [__p1B.get_child(0),__p1B.get_child(1)]
		__qCode = quarkPairs.get_quark_code(self.__maxPperpSquared,self.__activeQCodes) ##As maxPperpSquared now set to S123 of splitting.
		producedQuarkCodes.store(__qCode)
		__p1C0B, __p1C1B = __p1B.get_colour(0), __p1B.get_colour(1)
		__p1 = particles.particle(-__qCode,__p1B.get_four_momentum(),__p1BMs,__p1BCs,[0,0],__p1B.get_status_code())
		__p2 = particles.particle(__qCode,nV2,[0,0],[0,0],[0,0],1)
		__p3 = self.__chainList[self.__nextIndex][1].copy()
		__p1.set_four_momentum(nV1)
		__p3.set_four_momentum(nV3)
		__p1.set_unique_ID(pC.next()) ##Set unique particle IDs when produced.
		__p2.set_unique_ID(pC.next())
		__p1.set_produced_at(self.__maxPperpSquared) ##As now set to that of this produced particle.
		__p2.set_produced_at(self.__maxPperpSquared)
		__p1.set_mother(0,__p3.get_unique_ID())
		__p2.set_mother(1,__p3.get_unique_ID())
		if __p3.get_child(0) != 0: ##i.e don't overwrite.
			__p3.set_child(0,__p2.get_unique_ID())
		__p1, __p2, __p3 = self.set_new_colours_g_split(cC,__p1,__p2,__p3, __p1C0B, __p1C1B,"RHS")
		__nDPLHS = self.__chainList[self.__nextIndex].copy()
		__nDPRHS = self.__chainList[self.__nextIndex].copy()
		__nDPLHS[0], __nDPLHS[1] = __p1.copy(), __p2.copy()
		__nDPRHS[0], __nDPRHS[1] = __p2.copy(), __p3.copy()
		self.update_dipole_values(__p1.copy(),__nDPLHS,__nDPRHS,__p3.copy())
		self.__sideOfSplit = "RHS"

	def update_dipoles_g_split(self,pC,cC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a gluon is split."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		if (self.__dipoleRecoilIndex == 0): ##LHS recoiling.
			self.update_dipoles_g_split_LHS_recoil(pC,cC,nV1,nV2,nV3)
		else: ##RHS recoiling.
			self.update_dipoles_g_split_RHS_recoil(pC,cC,nV1,nV2,nV3)

	##~Photon Emission~##

	def update_dipole_values_p_prod(self,nPLHS,uDDP,nPRHS):
		"""A function to update the values in the surrounding dipoles for photon emission."""
		assert particles.check_is_particle(nPLHS)
		assert dipoles.check_is_dipole(uDDP)
		assert particles.check_is_particle(nPRHS)
		if ((not self.__isLoop) and (self.__nextIndex == 0)): ##No particle of dipole to left to update.
			pass
		else:
			__uLHIndex = self.get_chain_index(self.__nextIndex - 1)
			self.__chainList[__uLHIndex][1] = nPLHS.copy()
		if ((not self.__isLoop) and (self.__nextIndex == (len(self.__chainList) - 1))): ##-1 for list index from length.
			pass ##No particle of dipole to right to update.
		else:
			__uRHIndex = self.get_chain_index(self.__nextIndex + 1)
			self.__chainList[__uRHIndex][0] = nPRHS.copy()
		self.__chainList[self.__nextIndex] = uDDP.copy()

	def update_dipoles_p_prod_LHS_recoil(self,pC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a photon is emitted with the LHS recoiling."""
		assert counters.check_is_counter(pC)
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		__p1 = self.__chainList[self.__nextIndex][1].copy()
		__p2 = particles.particle('photon',nV2,[0,0],[0,0],[0,0],1)
		__p3 = self.__chainList[self.__nextIndex][0].copy()
		__p1.set_four_momentum(nV1)
		__p3.set_four_momentum(nV3)
		__p2.set_unique_ID(pC.next()) ##Set unique particle ID when produced.
		__p2.set_produced_at(self.__maxPperpSquared) ##As now set to that of this produced particle.
		__p2.set_mother(0,__p1.get_unique_ID())
		__p2.set_mother(1,__p3.get_unique_ID())
		if __p1.get_child(0) != 0: ##i.e don't overwrite
			__1.set_child(0,__p2.get_unique_ID())
		if __p3.get_child(1) != 0: ##i.e don't overwrite
			__p3.set_child(1,__p2.get_unique_ID())
		__uDDP = self.__chainList[self.__nextIndex].copy()
		__uDDP[0], __uDDP[1] = __p3.copy(), __p1.copy()
		self.update_dipole_values_p_prod(__p3.copy(),__uDDP,__p1.copy())
		self.__producedPhoton = __p2.copy()

	def update_dipoles_p_prod_RHS_recoil(self,pC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a photon is emitted with the RHS recoiling."""
		assert counters.check_is_counter(pC)
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		__p1 = self.__chainList[self.__nextIndex][0].copy()
		__p2 = particles.particle('photon',nV2,[0,0],[0,0],[0,0],1)
		__p3 = self.__chainList[self.__nextIndex][1].copy()
		__p1.set_four_momentum(nV1)
		__p3.set_four_momentum(nV3)
		__p2.set_unique_ID(pC.next()) ##Set unique particle ID when produced.
		__p2.set_produced_at(self.__maxPperpSquared) ##As now set to that of this produced particle.
		__p2.set_mother(0,__p1.get_unique_ID())
		__p2.set_mother(1,__p3.get_unique_ID())
		if __p1.get_child(1) != 0: ##i.e don't overwrite
			__p1.set_child(1,__p2.get_unique_ID())
		if __p3.get_child(0) != 0: ##i.e don't overwrite
			__p3.set_child(0,__p2.get_unique_ID())
		__uDDP = self.__chainList[self.__nextIndex].copy()
		__uDDP[0], __uDDP[1] = __p1.copy(), __p3.copy()
		self.update_dipole_values_p_prod(__p1.copy(),__uDDP,__p3.copy())
		self.__producedPhoton = __p2.copy()

	def update_dipoles_p_prod(self,pC,nV1,nV2,nV3):
		"""A function to update the relevant dipoles when a photon is emitted."""
		assert counters.check_is_counter(pC)
		for nV in [nV1,nV2,nV3]:
			assert fourVectors.check_is_fourVector(nV)
		assert (self.__maxPperpSquared > self.__cutOff)
		if (self.__dipoleRecoilIndex == 0): ##LHS recoiling.
			self.update_dipoles_p_prod_LHS_recoil(pC,nV1,nV2,nV3)
		else: ##RHS recoiling.
			self.update_dipoles_p_prod_RHS_recoil(pC,nV1,nV2,nV3)

	##~Control~##

	def evolve(self,activeQCodes,pC,cC):
		"""A function to evolve a chain through the next event if one is possible."""
		assert ((counters.check_is_counter(pC)) and (counters.check_is_counter(cC)))
		for __dipole in self.__chainList:
			assert __dipole.is_updated()
		assert (type(activeQCodes) == list)
		for qCode in activeQCodes:
			assert qCode in particleData.knownParticles.get_known_quarks()
		self.__activeQCodes = activeQCodes
		##If self.__maxPperpSquared > self.__cutOff, should get not self.__proceed.
		self.run_sudakovs() ##Get values and process for each dipole.
		self.__proceed = self.process_sudakovs() ##Work out what's happening next.
		if (not self.__proceed):
			self.set_showering_completed()
			return [0,None]
		else:
			self.prepare_dipoles()
			__nV1, __nV2, __nV3 = self.run_kinematics()
			if (self.__nextProcessCode == 1): ##Gluon emission.
				self.update_dipoles_g_prod(pC,cC,__nV1, __nV2, __nV3)
				return [1,None]
			elif (self.__nextProcessCode == 2): ##Gluon splitting.
				#Currently doesn't account for loops
				self.update_dipoles_g_split(pC,cC,__nV1, __nV2, __nV3)
				if (self.__sideOfSplit == "LHS"): ##i.e RHS is gluon splitting.
					__splitBeforeIndex = self.__nextIndex + 1
				elif (self.__sideOfSplit == "RHS"): ##i.e LHS is gluon splitting.
					__splitBeforeIndex = self.__nextIndex
				return [2,[__splitBeforeIndex,self.__nextPperpSquared]]
			elif (self.__nextProcessCode == 3): ##Photon emission.
				self.update_dipoles_p_prod(pC,__nV1,__nV2,__nV3)
				return [3,self.__producedPhoton]

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##
	import fourVectors
	import math

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "//////////////////////"
	print "Testing chains module:"
	print "//////////////////////"
	assertions.pause(__name__)
	
	##Setup here:##
	print "\nGenerating test values..."
	##Generate random test four-vectors:##
	testFourVectors = {'a':None,'b':None,'c':None,'d':None,'e':None}
	##Prevent the random vectors being the same as can't boost into frame of massless particle.
	testVector1, testVector2 = fourVectors.fourVector(0,0,0,0), fourVectors.fourVector(0,0,0,0) ##Required to start while loop.
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
	##Create test chains.
	testChain1 = chain([testParticle1,testParticle5])
	testChain2 = chain([testParticle1,testParticle2,testParticle3,testParticle4,testParticle5])
	##Also create loop test chain of centre gluons.
	testChain3 = chain([testParticle2,testParticle3,testParticle4],True)

	##Test check_is_chain:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_chain:\n"
	print "Calling check_is_chain on instance: " , check_is_chain(testChain2)
	print "Calling check_is_chain on second instance: " , check_is_chain(testChain3)
	print "Calling check_is_chain on wrong instance: " , check_is_chain(testVector1)
	print "Calling check_is_chain on second wrong instance: " , check_is_chain(testChain3[0])
	print "Calling check_is_chain on 1.055: " , check_is_chain(1.055)
	print "Calling check_is_chain on 'word': " , check_is_chain('word')
	results1 = [check_is_chain(testChain2),check_is_chain(testChain3),check_is_chain(testVector1)]
	results2 = [check_is_chain(testChain3[0]),check_is_chain(1.055),check_is_chain('word')]
	results = results1 + results2
	if ((sum(results) == 2) and (results[0] == True) and (results[1] == True)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_chain."
	assertions.pause(__name__)

	##Test chain class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////"
	print "Testing chain class:"
	print "////////////////////"
	assertions.pause(__name__)

	##__init__() tested implicitly.

	##Test __str__(), simple_str() and __repr__() functions:##
	print "\n--------------------------------------------------\n"
	print "Testing __str__(), simple_str() and __repr__() functions:\n"
	for testChain in [testChain1,testChain2,testChain3]:
		print "String:\n", testChain
		assertions.pause(__name__)
		print "\nReprsentation:\n", repr(testChain)
		assertions.pause(__name__)
		print "\nSimple string:", testChain.simple_str()
		print "-----"
		assertions.pause(__name__)
	print "\nFinished testing __str__(), simple_str() and __repr__() functions."
	assertions.pause(__name__)

	##Test check_is_closed function:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_closed function:\n"
	print "Calling check is closed on the three chains above returns:"
	for testChain in [testChain1,testChain2,testChain3]:
		print testChain.check_is_closed()
	if ((not testChain1.check_is_closed()) and (not testChain2.check_is_closed()) and testChain3.check_is_closed()):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_closed function."
	assertions.pause(__name__)

	##Test get_chain_index, get_chain_list and __getitem__ functions:##
	print "\n--------------------------------------------------\n"
	print "Testing get_chain_index, get_chain_list and __getitem__ functions:\n"
	for testChain in [testChain1,testChain2,testChain3]:
		counter = 1
		print "\nUsing the chain:\n", testChain.simple_str(), "\n"
		for anIndex in range(-10,11):
			print "Calling with index", anIndex, "returns:", testChain.get_chain_index(anIndex)
			print "which gives the dipole:", testChain[anIndex].simple_str()
			if counter%3 == 0:
				assertions.pause(__name__)
			counter += 1
		print "\nFinished checking for", testChain.simple_str()
		assertions.pause(__name__)
		print "\nget_chain_list returns:\n\n", testChain.get_chain_list()
		assertions.pause(__name__)
	print "\nFinished testing get_chain_index, get_chain_list and __getitem__ functions."
	assertions.pause(__name__)

	##Test get_max_Pperp_squared and get_next_Pperp_squared functions:##
	print "\n--------------------------------------------------\n"
	print "Testing get_max_Pperp_squared and get_next_Pperp_squared functions:\n"
	for testChain in [testChain1,testChain2,testChain3]:
		print "Using", testChain.simple_str()
		print "Calling get_max_Pperp_squared returns:", testChain.get_max_Pperp_squared()
		print "Calling get_next_Pperp_squared returns:", testChain.get_next_Pperp_squared()
	print "\nFinished testing get_max_Pperp_squared and get_next_Pperp_squared functions."
	assertions.pause(__name__)

	##Test evolve() and all sub-functions for gluon emission:##
	print "\n--------------------------------------------------\n"
	print "Testing evolve() and all sub-functions for gluon emission:\n"
	for testChain in [testChain1,testChain2,testChain3]:
		print "\nUsing the chain:\n", testChain.simple_str(), "\n"
		for n in range(20):
			print "Has maxPperpSquared:", testChain.get_max_Pperp_squared()
			testChain.evolve()
			print "\nHad next PperpSquared:", testChain.get_next_Pperp_squared()
			print "\n", testChain.simple_str(), "\n"
			print str(testChain)
			assertions.pause(__name__)
		print "\nFinished checking for", testChain.simple_str()
		assertions.pause(__name__)
	print "\nFinished testing evolve() and all sub-functions for gluon emission."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "////////////////////////////////"
	print "Finished checking chains module!"
	print "////////////////////////////////"