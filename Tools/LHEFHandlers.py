####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module for reading and writing LHEF XML files."""

##Uses DOM instead of SAX which was the other option.
##Document Object Model should be faster but could use more memory.

##Import required modules:##
import xml.dom.minidom
import os
import pyperclip
import datetime
import assertions
import precision
import constants
import particleData
import fourVectors
import particles
import kinematics

print "\n////////////////////////////"
print "Loading LHEFHandlers module:"
print "Test code not written yet!"
print "Init values not coded in yet!"
print "////////////////////////////\n"
assertions.pause_loading_module() #can move import lower when removed here.

##Functions:##

def check_is_LHEFReader(toCheck):
	"""A function to check for an instance of the LHEFReader class."""
	return isinstance(toCheck,LHEFReader)

def check_is_LHEFShowerWriter(toCheck):
	"""A function to check for an instance of the LHEFShowerWriter class."""
	return isinstance(toCheck,LHEFShowerWriter)

def check_is_LHEFMEWriter(toCheck):
	"""A function to check for an instance of the LHEFMEWriter class."""
	return isinstance(toCheck,LHEFMEWriter)

def zero_or(value):
	"""A function to see if a value can be set to 0.0 when writing an LHEF file."""
	##Main use is to remove -ve energies that are actually zero with the code precision.
	assert (type(value) == float)
	if precision.check_numbers_equal(value,0.0):
		return 0.0
	else:
		return value

##Classes:##

class LHEFReader(object):
	"""A class for reading events from LHEF XML files."""

	def __init__(self,fileName):
		"""A function to initiate an LHEF XML file reader."""
		assert (type(fileName) == str)
		if not (fileName[-4:] == '.lhe'):
			self.__toRead = fileName + '.lhe'
		else:
			self.__toRead = fileName
		##Open XML document using minidom parser.
		print "\nReading file...\n"
		self.__DOMTree = xml.dom.minidom.parse(self.__toRead)
		self.__data = self.__DOMTree.documentElement
		self.__events = self.__data.getElementsByTagName("event")

	def get_source_string(self):
		"""A function to return a string of the source file used."""
		return self.__toRead

	def get_all_as_string(self):
		"""A function to return the whole LHEF XML file as a string."""
		with open(self.__toRead, 'r') as fileToRead:
			__theString = fileToRead.read() #.replace('\n', '')
		return __theString

	def get_number_events(self):
		"""A function to return the number of events in a LHEF XML file."""
		return len(self.__events)

	def get_event_object(self,eventIndex):
		"""A function to return an event object from a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.__events[eventIndex]

	def get_event_trials(self,eventIndex):
		"""A function to get the number of trials for an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		__theEventObj = self.get_event_object(eventIndex)
		return __theEventObj.getAttribute("trials")

	def get_event_muf2(self,eventIndex):
		"""A function to get the muf2 value for an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		__theEventObj = self.get_event_object(eventIndex)
		return __theEventObj.getAttribute("muf2")

	def get_event_mur2(self,eventIndex):
		"""A function to get the mur2 value for an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		__theEventObj = self.get_event_object(eventIndex)
		return __theEventObj.getAttribute("mur2")

	def get_event_row_strings(self,eventIndex):
		"""A function to return the data of an event from a LHEF XML file."""
		##Returns a list of each line in the event, ignoring blank lines.
		assert (type(eventIndex) == int)
		__theEventStr = str(self.__events[eventIndex].childNodes[0].data)
		__theEventLines = __theEventStr.split('\n')
		__theEventLines = __theEventLines[1:-1] ##Cut off blank lines.
		return __theEventLines

	def get_event_top_row_values(self,eventIndex):
		"""A function to return the top row of a LHEF XML file."""
		assert (type(eventIndex) == int)
		__topRowValues = self.get_event_row_strings(eventIndex)[0]
		##Strip out any spaces
		__topRowValues = __topRowValues.replace('\t',' ').split(' ')
		return [x for x in __topRowValues if (x != '' and x != '\t')]

	def get_event_number_particles(self,eventIndex):
		"""A function to get the number of particles in an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.get_event_top_row_values(eventIndex)[0]

	def get_event_process_ID(self,eventIndex):
		"""A function to get the process ID in an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.get_event_top_row_values(eventIndex)[1]

	def get_event_weight(self,eventIndex):
		"""A function to get the weight in an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.get_event_top_row_values(eventIndex)[2]

	def get_event_scale(self,eventIndex):
		"""A function to get the scale in an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.get_event_top_row_values(eventIndex)[3]

	def get_event_alphaEM(self,eventIndex):
		"""A function to get the EM coupling constant in an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.get_event_top_row_values(eventIndex)[4]

	def get_event_alphaS(self,eventIndex):
		"""A function to get the strong coupling constant in an event in a LHEF XML file."""
		assert (type(eventIndex) == int)
		return self.get_event_top_row_values(eventIndex)[5]

	def get_event_particle_strings(self,eventIndex,particleIndex):
		"""A function to return the details of a particle from an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		##Ignore the top line which is not a particle gives +1 in index below.
		__theDetails = self.get_event_row_strings(eventIndex)[particleIndex + 1].replace('\t',' ').split(' ')
		__theDetails = [x for x in __theDetails if (x != '' and x != '\t')]
		return __theDetails

	def get_event_particle_ID(self,eventIndex,particleIndex):
		"""A function to get the ID of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		return int(self.get_event_particle_strings(eventIndex,particleIndex)[0])

	def get_event_particle_when(self,eventIndex,particleIndex):
		"""A function to get the initial/intermediate/final state code of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		##Assumed when code should be an integer when outputting.
		return int(self.get_event_particle_strings(eventIndex,particleIndex)[1])

	def get_event_particle_mothers(self,eventIndex,particleIndex):
		"""A function to get the mothers of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		##Assumed mother codes should be an integer when outputting.
		__mother1 = int(self.get_event_particle_strings(eventIndex,particleIndex)[2])
		__mother2 = int(self.get_event_particle_strings(eventIndex,particleIndex)[3])
		return [__mother1,__mother2]

	def get_event_particle_colours(self,eventIndex,particleIndex):
		"""A function to get the coloursL of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		##Assumed colour codes should be an integer when outputting.
		__colour1 = int(self.get_event_particle_strings(eventIndex,particleIndex)[4])
		__colour2 = int(self.get_event_particle_strings(eventIndex,particleIndex)[5])
		return [__colour1,__colour2]

	def get_event_particle_four_momentum(self,eventIndex,particleIndex):
		"""A function to get the four-momentum of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		__fourVector = fourVectors.fourVector()
		##Energy is the last component in the LHEF
		__fourVector[0] = float(self.get_event_particle_strings(eventIndex,particleIndex)[9])
		__fourVector[1] = float(self.get_event_particle_strings(eventIndex,particleIndex)[6])
		__fourVector[2] = float(self.get_event_particle_strings(eventIndex,particleIndex)[7])
		__fourVector[3] = float(self.get_event_particle_strings(eventIndex,particleIndex)[8])
		return __fourVector

	def get_event_particle_mass(self,eventIndex,particleIndex):
		"""A function to get the mass of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		return float(self.get_event_particle_strings(eventIndex,particleIndex)[10])

	def get_event_particle_lifetime(self,eventIndex,particleIndex):
		"""A function to get the proper lifetime of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		return float(self.get_event_particle_strings(eventIndex,particleIndex)[11])

	def get_event_particle_spin(self,eventIndex,particleIndex):
		"""A function to get the spin of a particle in an event in an LHEF XML file."""
		assert ((type(eventIndex) == int) and (type(particleIndex) == int))
		return float(self.get_event_particle_strings(eventIndex,particleIndex)[12])

	def get_event_particle(self,eventIndex,particleIndex):
		"""A function to return a particle object from an LHEF XML file given it's indices."""
		__code = self.get_event_particle_ID(eventIndex,particleIndex)
		__fM = self.get_event_particle_four_momentum(eventIndex,particleIndex)
		__mothers = self.get_event_particle_mothers(eventIndex,particleIndex)
		__children = [0,0] ##LHEF doesn't have children; only mothers.
		__colours = self.get_event_particle_colours(eventIndex,particleIndex)
		__statusCode = self.get_event_particle_when(eventIndex,particleIndex)
		__theParticle = particles.particle(__code,__fM,__mothers,__children,__colours,__statusCode)
		##Assuming the LHEF file is appropriately ordered:
		__theParticle.set_unique_ID(particleIndex+1) ##+1 as start at particle 1.
		return __theParticle

class LHEFShowerWriter(object):
	"""A class for writing LHEF XML files for Dipole Showered events."""

	def __init__(self,sourceFileReader,numberEvents):
		"""A function to initialise a Dipole Shower LHEF XML writer."""
		assert check_is_LHEFReader(sourceFileReader)
		assert (type(numberEvents) == int)
		self.__numberEvents = numberEvents
		##Strips out and replaces events in source LHEF file, only changing required values.
		self.__readFrom = sourceFileReader
		self.__sourceFileString = self.__readFrom.get_source_string()
		self.__stringFile = self.__readFrom.get_all_as_string()
		self.__top = self.__stringFile.split('\n<event')[0]
		self.__bottom = "\n</LesHouchesEvents>"
		self.__numberEventsSaved = 0
		self.__haveEvents = False ##To track if at least one has been added.
		self.__topWritten = False

	def write_top(self):
		"""A function to update the sections of a LHEF file above the events."""
		if not "#Showered by PyShower" in self.__top:
			__splitTop = self.__top.split("-->")
			self.__top = __splitTop[0] + "\n##Showered by PyShower (Daniel Osborne 2015/16 - Durham University)\n"
			self.__top += "##ME LHEF source file: " + self.__sourceFileString + "\n\n-->" + __splitTop[1]
		self.create()
		self.__topWritten = True

	def add_event(self,trials,muf2,mur2,processID,weight,scale,alphaEM,alphaS,particlesIn,particlesOut):
		"""A function to convert a produced particles data to LHEF and add to the writer."""
		assert ((type(trials) == int) or (isinstance(trials, basestring)))
		assert ((type(muf2) == int) or (type(muf2) == float) or (isinstance(muf2, basestring)))
		assert ((type(mur2) == int) or (type(mur2) == float) or (isinstance(mur2, basestring)))
		assert ((type(processID) == int) or (isinstance(processID, basestring)))
		assert ((type(weight) == float) or (type(weight) == int) or (isinstance(weight, basestring)))
		assert ((type(scale) == float) or (isinstance(scale, basestring)))
		assert ((type(alphaEM) == float) or (isinstance(alphaEM, basestring)))
		assert ((type(alphaS) == float) or (isinstance(alphaS, basestring)))
		if (self.__topWritten == False):
			self.write_top()
			self.__haveEvents = True ##Only called once here.
		__listOfParticles = particlesIn + particlesOut
		#Currently not sorting by order of production.
		#__listOfIDs = [x.get_unique_ID() for x in __listOfParticles]
		#__sortedListOfParticles = [x for (y,x) in sorted(zip(__listOfIDs,__listOfParticles), key=lambda pair: pair[0])]
		__sortedListOfParticles = __listOfParticles ##Replacement line to keep the chain ordering.
		for particle in __sortedListOfParticles:
			assert particles.check_is_particle(particle)
			assert particle.__nonzero__()
		__numberParticles = len(particlesIn) + len(particlesOut)
		__line1 = "\n<event>" # trials='" + str(trials) + "' muf2='" + str(muf2) + "' mur2='" + str(mur2) + "'>"
		__line2 = ("\n\t\t" + str(__numberParticles) + "\t" + str(processID) + "\t" + str(weight) + "\t" + str(scale))
		__line2 += ("\t" + str(alphaEM) + "\t" + str(alphaS))
		__inLines, __outLines = "", ""
		__endLine = "\n</event>"
		__zO = zero_or
		__pHS = "{:." + str(constants.LHEF_DS_number_decimal_places()) + "e}"
		for i, p in enumerate(__sortedListOfParticles): ##Iterate through both at once, but with particlesIn first!
			if not precision.check_numbers_equal(p.get_four_momentum()*p.get_four_momentum(),0.0):
				fFV = kinematics.fix_not_massless(p.get_four_momentum()) ##Adjust to make sure they are exactly massless for Pythia.
			else:
				fFV = p.get_four_momentum()
			if i <= 1: ##Carry on as normal for first two (initial-state) particles.
				__mum1, __mum2 = p.get_mother(0), p.get_mother(1)
			else: ##Set all mother values to 1 and 2 as required by Pythia 8.2.
				__mum1, __mum2 = 1, 2
			__newLine = ("\n\t\t\t" + str(p.get_code()) + "\t" + str(p.get_status_code()) + "\t" + str(__mum1))
			__newLine += ("\t" + str(__mum2) + "\t" + str(p.get_colour(0)) + "\t\t" + str(p.get_colour(1)))
			__newLine += ("\t\t" + __pHS.format(__zO(fFV[1])) + "\t" + __pHS.format(__zO(fFV[2])) + "\t")
			__newLine += (__pHS.format(__zO(fFV[3])) + "\t" + __pHS.format(__zO(fFV[0])) + "\t")
			__newLine += (__pHS.format(__zO(p.get_mass())) + "\t" + str(__zO(p.get_width())) + "\t\t" + str(p.get_spin()))
			__inLines += __newLine
		__toSave = __line1 + __line2 + __inLines + __outLines + __endLine
		__file = open(self.__saveName,'a') ##'a' means append to.
		__file.write(__toSave)
		__file.close()
		self.__numberEventsSaved += 1

	def create(self):
		"""A function to create the LHEF XML file."""
		assert not self.__topWritten
		##Get today's date for results folder name.
		__date = datetime.datetime.now().strftime("%d.%m.%y")
		if ("MEs_" in self.__sourceFileString):
			__split1 = self.__sourceFileString.split('MEs_')
			__split2 = __split1[0].split('\\')
			__saveName = __split2[2] + "_DS" + datetime.datetime.now().strftime("at%H.%M.%S") + '.lhe'
		else:
			__split3 = self.__sourceFileString.split('\\')
			__saveName = 'DS' + datetime.datetime.now().strftime("at%H.%M.%S") + "_" + __split3[2]
		if not (__saveName[-4:] == '.lhe'):
			__saveName = __saveName + '.lhe'
		__saveName = "ME_DSs\\" + __date + "\\" + __saveName
		__toSave = self.__top
		##Make the directory for the date if it doesn't exist.
		if not os.path.exists(os.path.dirname(__saveName)):
			try:
				os.makedirs(os.path.dirname(__saveName))
			except OSError as exc: ##Guard against race condition.
				if exc.errno != errno.EEXIST:
					raise
		__file = open(__saveName,'w') ##'w' will create new file but also write over if already existing!
		__file.write(__toSave)
		__file.close()
		self.__saveName = __saveName ##So it can be used elsewhere to access the file.

	def save(self):
		"""A function to finalise the produced LHEF XML file."""
		assert self.__haveEvents
		if not (self.__numberEvents == self.__numberEventsSaved):
			print "\nWARNING: Missing" + str( self.__numberEvents - self.__numberEventsSaved) + "events.\n"
		__toSave = self.__bottom
		__file = open(self.__saveName,'a') ##'a' means append to.
		__file.write(__toSave)
		__file.close()
		print "\n////////////////////////////////////////////////////////////////////"
		print "Saved as:", self.__saveName
		print "////////////////////////////////////////////////////////////////////\n"
		pyperclip.copy(self.__saveName)
		print "\nThis file path has been copied to your clipboard for ease of use."

class LHEFMEWriter(object):
	"""A class for writing LHEF XML files for ME events."""

	def __init__(self,processName,numberEvents,S123,qCodes):
		"""A function to initialise a ME LHEF XML writer."""
		knownQuarkCodes = particleData.knownParticles.get_known_quarks()
		assert (type(processName) == str)
		assert (type(numberEvents) == int)
		assert ((type(S123) == float) or (isinstance(S123, basestring)))
		assert (type(qCodes) == list)
		for qCode in qCodes:
			assert ((type(qCode) == int) and (qCode in knownQuarkCodes))
		self.__numberEvents, self.__S123, self.__qCodes = numberEvents, S123, qCodes
		##Strips out and replaces events in source LHEF file, only changing required values.
		self.__processName = processName
		__templateLoc = "Templates/LHEF_ME_Template.lhe"
		self.__templateFile = LHEFReader(__templateLoc)
		self.__stringFile = self.__templateFile.get_all_as_string()
		self.__top = self.__stringFile.split('\n<event')[0]
		self.__bottom = "\n</LesHouchesEvents>"
		self.__numberEventsSaved = 0
		self.__haveEvents = False ##To track if at least one has been added.
		self.__topWritten = False

	def write_top(self):
		"""A function to update the sections of a LHEF file above the events."""
		__splitTop = self.__top.split("-->")
		self.__top = __splitTop[0]
		if not "#Matrix elements approximated by PyShower" in self.__top:
			self.__top += ("\n##Matrix elements approximated by PyShower (Daniel Osborne 2015/16 - Durham University)\n")
		self.__top += ("##numEvents='" + str(self.__numberEvents) + "', S123='" + str(self.__S123) + "', qCodes=" + str(self.__qCodes))
		self.__top += "\n\n-->" + __splitTop[1]
		self.create()
		self.__topWritten = True

	def add_event(self,trials,muf2,mur2,processID,weight,scale,alphaEM,alphaS,particlesIn,particlesOut):
		"""A function to convert a produced particles data to LHEF and add to the writer."""
		assert ((type(trials) == int) or (isinstance(trials, basestring)))
		assert ((type(muf2) == int) or (type(muf2) == float) or (isinstance(muf2, basestring)))
		assert ((type(mur2) == int) or (type(mur2) == float) or (isinstance(mur2, basestring)))
		assert ((type(processID) == int) or (isinstance(processID, basestring)))
		assert ((type(weight) == float) or (type(weight) == int) or (isinstance(weight, basestring)))
		assert ((type(scale) == float) or (isinstance(scale, basestring)))
		assert ((type(alphaEM) == float) or (isinstance(alphaEM, basestring)))
		assert ((type(alphaS) == float) or (isinstance(alphaS, basestring)))
		if (self.__topWritten == False):
			self.write_top()
			self.__haveEvents = True ##Only called once here.
		__listOfParticles = particlesIn + particlesOut
		__listOfIDs = [x.get_unique_ID() for x in __listOfParticles]
		__sortedListOfParticles = [x for (y,x) in sorted(zip(__listOfIDs,__listOfParticles), key=lambda pair: pair[0])]
		for particle in __sortedListOfParticles:
			assert particles.check_is_particle(particle)
			assert particle.__nonzero__()
		__numberParticles = len(particlesIn) + len(particlesOut)
		__line1 = "\n<event>" #trials='" + str(trials) + "' muf2='" + str(muf2) + "' mur2='" + str(mur2) + "'>"
		__line2 = ("\n\t\t" + str(__numberParticles) + "\t" + str(processID) + "\t" + str(weight) + "\t" + str(scale))
		__line2 += ("\t" + str(alphaEM) + "\t" + str(alphaS))
		__inLines, __outLines = "", ""
		__endLine = "\n</event>"
		__zO = zero_or
		__pHS = "{:." + str(constants.LHEF_ME_number_decimal_places()) + "e}"
		for p in __sortedListOfParticles: ##Iterate through both at once, but with particlesIn first!
			if not precision.check_numbers_equal(p.get_four_momentum()*p.get_four_momentum(),0.0): ##Shouldn't ever need to call it.
				fFV = kinematics.fix_not_massless(p.get_four_momentum()) ##Adjust to make sure they are exactly massless for Pythia.
			else:
				fFV = p.get_four_momentum()
			__newLine = ("\n\t\t\t" + str(p.get_code()) + "\t" + str(p.get_status_code()) + "\t" + str(p.get_mother(0)))
			__newLine += ("\t" + str(p.get_mother(1)) + "\t" + str(p.get_colour(0)) + "\t\t" + str(p.get_colour(1)))
			__newLine += ("\t\t" + __pHS.format(__zO(fFV[1])) + "\t" + __pHS.format(__zO(fFV[2])) + "\t")
			__newLine += (__pHS.format(__zO(fFV[3])) + "\t" + __pHS.format(__zO(fFV[0])) + "\t")
			__newLine += (__pHS.format(__zO(p.get_mass())) + "\t" + str(__zO(p.get_width())) + "\t\t" + str(p.get_spin()))
			__inLines += __newLine
		__toSave = __line1 + __line2 + __inLines + __outLines + __endLine
		__file = open(self.__saveName,'a') ##'a' means append to.
		__file.write(__toSave)
		__file.close()
		self.__numberEventsSaved += 1

	def create(self):
		"""A function to create the LHEF XML file."""
		assert not self.__topWritten
		##Get today's date for results folder name.
		__date = datetime.datetime.now().strftime("%d.%m.%y")
		__saveName = self.__processName + "_" + str(self.__numberEvents) + 'MEs_' + datetime.datetime.now().strftime("at%H.%M.%S") + '.lhe'
		__saveName = "MEs\\" + __date + "\\" + __saveName
		__toSave = self.__top
		##Make the directory for the date if it doesn't exist.
		if not os.path.exists(os.path.dirname(__saveName)):
			try:
				os.makedirs(os.path.dirname(__saveName))
			except OSError as exc: ##Guard against race condition.
				if exc.errno != errno.EEXIST:
					raise
		__file = open(__saveName,'w') ##'w' will create new file but also write over if already existing!
		__file.write(__toSave)
		__file.close()
		self.__saveName = __saveName ##So it can be used elsewhere to access the file.

	def save(self):
		"""A function to finalise the produced LHEF XML file."""
		assert self.__haveEvents
		if not (self.__numberEvents == self.__numberEventsSaved):
			print "\nWARNING: Missing" + str(self.__numberEvents - self.__numberEventsSaved) + "events.\n"
		__toSave = self.__bottom
		__file = open(self.__saveName,'a') ##'a' means append to.
		__file.write(__toSave)
		__file.close()
		print "\n////////////////////////////////////////////////////////////////////"
		print "Saved as:", self.__saveName
		print "////////////////////////////////////////////////////////////////////\n"
		pyperclip.copy(self.__saveName)
		print "\nThis file path has been copied to your clipboard for ease of use."

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "////////////////////////////"
	print "Testing LHEFHandlers module:"
	print "////////////////////////////"
	assertions.pause(__name__)

	##Setup here:##
	testPathIn = "../MEs/belle.lhe"
	testReader = LHEFReader(testPathIn)
	testPathOut = "testLHEFOut.lhe"
	testWriter = LHEFShowerWriter(testPathIn)

	##Test check_is_LHEFReader and check_is_LHEFShowerWriter:##
	print "\n--------------------------------------------------\n"
	print "Testing check_is_LHEFReader and check_is_LHEFShowerWriter:"
	print "\nCalling check_is_LHEFReader on instance: " , check_is_LHEFReader(testReader)
	print "Calling check_is_LHEFReader on wrong instance: " , check_is_LHEFReader(testWriter)
	print "Calling check_is_LHEFReader on 1.055: " , check_is_LHEFReader(1.055)
	print "Calling check_is_LHEFReader on 'word': " , check_is_LHEFReader('word')
	testResults = [check_is_LHEFReader(testReader),check_is_LHEFReader(testWriter)]
	testResults += [check_is_LHEFReader(1.055),check_is_LHEFReader('word')]
	print "Calling check_is_LHEFShowerWriter on instance: " , check_is_LHEFShowerWriter(testWriter)
	print "Calling check_is_LHEFShowerWriter on wrong instance: " , check_is_LHEFShowerWriter(testReader)
	print "Calling check_is_LHEFShowerWriter on 1.055: " , check_is_LHEFShowerWriter(1.055)
	print "Calling check_is_LHEFShowerWriter on 'word': " , check_is_LHEFShowerWriter('word')
	testResults += [check_is_LHEFShowerWriter(testWriter),check_is_LHEFShowerWriter(testReader)]
	testResults += [check_is_LHEFShowerWriter(1.055),check_is_LHEFShowerWriter('word')]
	if (sum(testResults) == 2):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_is_LHEFReader and check_is_LHEFShowerWriter."
	assertions.pause(__name__)

	##Test LHEFReader class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////////"
	print "Testing LHEFReader class:"
	print "/////////////////////////"
	assertions.pause(__name__)

	##__init__() tested implicitly.

	##Test get_all_as_string and get_number_events functions:##
	print "\n--------------------------------------------------\n"
	print "Testing get_all_as_string and get_number_events functions:\n"
	print "Calling for " + testPathIn + ":"
	assertions.pause(__name__)
	print testReader.get_all_as_string()
	assertions.pause(__name__)
	print "Calling get_number_events returns:", testReader.get_number_events()
	if (testReader.get_number_events() == 100):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing get_all_as_string and get_number_events functions."
	assertions.pause(__name__)

	##get_event_object tested implicitly below:
	##Test get_event_trials, get_event muf2, get_event_mur2 functions:##
	print "\n--------------------------------------------------\n"
	print "Testing get_event_trials, get_event muf2, get_event_mur2 functions:\n"
	print "Using the first event in " + testPathIn + ":"
	print "Calling get_event_trials returns:", testReader.get_event_trials(0)
	print "Calling get_event_muf2 returns:", testReader.get_event_muf2(0)
	print "Calling get_event_mur2 returns:", testReader.get_event_mur2(0)
	results = [testReader.get_event_trials(0),testReader.get_event_muf2(0),testReader.get_event_mur2(0)]
	expected = [2,8317.44,8317.44]
	if (results == expected):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing get_event_trials, get_event muf2, get_event_mur2 functions."
	assertions.pause(__name__)

	##get_event_row_strings and get_event_top_row_values tested implicitly below
	##Test get_event _number_particles, _process_ID, _weight, _scale, _alphaEM and _alphaS functions:##
	print "\n--------------------------------------------------\n"
	print "Testing get_event _number_particles, _process_ID, _weight, _scale, _alphaEM and _alphaS functions:\n"
	print "Using the first event in " + testPathIn + ":"
	print "Calling get_event_number_particles returns:", testReader.get_event_number_particles(0)
	print "Calling get_event_process_ID returns:", testReader.get_event_process_ID(0)
	print "Calling get_event_weight returns:", testReader.get_event_weight(0)
	print "Calling get_event_scale returns:", testReader.get_event_scale(0)
	print "Calling get_event_alphaEM returns:", testReader.get_event_alphaEM(0)
	print "Calling get_event_alphaS returns:", testReader.get_event_alphaS(0)
	results = [testReader.get_event_number_particles(0),testReader.get_event_process_ID(0),testReader.get_event_weight(0)]
	results += [testReader.get_event_scale(0),testReader.get_event_alphaEM(0),testReader.get_event_alphaS(0)]
	expected = [4,1,1.2798299414e+04,9.1200000000e+01,-1.0000000000e+00,1.1879754368e-01]
	if (results == expected):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing get_event _number_particles, _process_ID, _weight, _scale, _alphaEM and _alphaS functions."
	assertions.pause(__name__)

	##Test all get_event_particle_'' functions:##
	print "\n--------------------------------------------------\n"
	print "Testing all get_event_particle_'' functions:\n"
	print "Using the first event in '" + testPathIn + "':"
	print "\nThe PDGID is:", testReader.get_event_particle_ID(0,0)
	print "The 'when' code is:", testReader.get_event_particle_when(0,0)
	print "The mothers are:", testReader.get_event_particle_mothers(0,0)
	print "The coloursL are:", testReader.get_event_particle_colours(0,0)
	print "The four-momentum is:", testReader.get_event_particle_four_momentum(0,0)
	print "The mass is:", testReader.get_event_particle_mass(0,0)
	print "The proper lifetime is:", testReader.get_event_particle_lifetime(0,0)
	print "The spin is:", testReader.get_event_particle_spin(0,0)
	results = [testReader.get_event_particle_ID(0,0),testReader.get_event_particle_when(0,0),testReader.get_event_particle_mothers(0,0)]
	results += [testReader.get_event_particle_colours(0,0),testReader.get_event_particle_four_momentum(0,0),testReader.get_event_particle_mass(0,0)]
	results += [testReader.get_event_particle_lifetime(0,0),testReader.gget_event_particle_spin(0,0)]
	expected = [11,-1,0,0,0,0,fourVectors.fourVector(5.2899734145e+00,0.0000000000e+00,0.0000000000e+00,5.2899734145e+00),0,0,9]
	if (results == expected):
		print "Test successful!"
	else:
		print "Test failed!"
	print "\nFinished testing all get_event_particle_'' functions."
	assertions.pause(__name__)

	##Test LHEFShowerWriter class:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "///////////////////////////////"
	print "Testing LHEFShowerWriter class:"
	print "///////////////////////////////"
	assertions.pause(__name__)

#	##Test __str__(), simple_str() and __repr__() functions:##
#	print "\n--------------------------------------------------\n"
#	print "Testing __str__(), simple_str() and __repr__() functions:\n"
#	print "Calling in order on the three test results containers used above returns:"
#	for testRC in testRCs:
#		print "\n", str(testRC)
#		print "\n", testRC.simple_str()
#		print "\n", repr(testRC)
#		assertions.pause(__name__)
#	print "\nFinished testing __str__(), simple_str() and __repr__() functions."
#	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "//////////////////////////////////////"
	print "Finished checking LHEFHandlers module!"
	print "//////////////////////////////////////"