####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module containing functions for carrying out checks during dipole showering."""

##Import required modules:##
import numpy
from matplotlib import pyplot

print "\n//////////////////////////"
print "Loading assertions module:"
print "//////////////////////////\n"

##Functions:##

##~~Control Values~~#

def gluon_splitting_on():
	"""A function to set whether gluon splitting to qqBar pairs is turned on."""
	##Controlled by changing which line below is/isn't commented out.
	#return False
	return True

def EM_radiation_on():
	"""A function to set whether photon emission from qqBar pairs is turned on."""
	##Controlled by changing which line below is/isn't commented out.
	#return False
	return True

def pause(currentName):
	"""A function to control pausing the running code during testing."""
	##Swap the results to determine wether to run with or without pausing.
	if currentName == "__main__": ##A pause in test code.
		raw_input("\nHit enter to continue: ")
		#pass
	else: ##In main code.
		raw_input("\nHit enter to continue: ")
		#pass

def pause_loading_module():
	"""A function to control pausing a loading module if it has important reminders."""
	##Swap the results to determine wether to run with or without pausing.
	raw_input("\nHit enter to continue: ")
	#pass

def show_graph():
	"""A function to control displaying graphs during testing."""
	##Swap the results to determine wether graphs show up or not.
	pyplot.show()
	#pass

##~~End Control Values~~#

def all_are_numbers(valuesToCheck):
	"""A function to check whether all the values in the list/tuple passed are numbers of the correct type."""
	__typesCorrect = True
	for v in valuesToCheck:
		if (type(v) != float):
			if (type(v) != int):
				if (type(v) != long):
					if (type(v) != numpy.float64):
						__typesCorrect = False
	return __typesCorrect

def force_float_number(number):
	"""A function to float an integer but leave doubles."""
	##Used to prevent incorrect division by float errorsself.
	assert all_are_numbers([number])
	return float(number)

def check_float(potentialFloat):
	"""A function to check if an object is a float."""
	return (type(potentialFloat) == float)

def valid_four_vector_index(index):
	"""A function to check an index is within the fourVector range: 0 <= i <= 3."""
	if (type(index) == int):
		if ((0 <= index) and (index <= 3)):
			return True
	return False

def valid_dipole_index(index):
	"""A function to check an index is within the dipole range: 0 <= i <= 1."""
	if (type(index) == int):
		if ((0 <= index) and (index <= 1)):
			return True
	return False

def valid_particle_status_code(index):
	"""A function to check an index is a valid particle status code."""
	if index in [-1,0,1]:
		return True
	else:
		return False

def safe_division(denominator):
	"""A function to check that a denominator is not zero within the code precision."""
	assert all_are_numbers([denominator])
	if precision.division_check_numbers_equal(denominator,0.0):
		return False
	else:
		return True

##Classes:##

##Import here to prevent cyclic failure:
import precision

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##

	##Begin testing:##
	print " \n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "//////////////////////////"
	print "Testing assertions module:"
	print "//////////////////////////"
	pause(__name__)

	##Setup here:##
	print "\nGenerating test values..."
	tAllNumbers = [12,20,57,20.025,6.0202,52585.,245435,3535.668868,64684,float(1)]
	tNotAllNumbers = [12,20,'word',20.025,6.0202,[52585],long(245435),long(3535.668868),64684,float(1)]
	tAnInteger = 238
	tAFloat = float(253.658)
	tALong = long(586.359)
	tIndices = [-2,-1,0,1,2,3,4,5,6,"bin",658,2.0]
	tDenominators = [1,2,86.0,6.664,0.000055,0,65,0.0000,-254,-3.564,6.37*(10**(-12)),6.37*(10**(-10))]

	##Test gluon_splitting_on:##
	print "\n--------------------------------------------------\n"
	print "Testing gluon_splitting_on:\n"
	print "Calling it returns: " , gluon_splitting_on()
	if (gluon_splitting_on() or (not gluon_splitting_on())):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing gluon_splitting_on."
	pause(__name__)

	##Test EM_radiation_on:##
	print "\n--------------------------------------------------\n"
	print "Testing EM_radiation_on:\n"
	print "Calling it returns: " , EM_radiation_on()
	if (EM_radiation_on() or (not EM_radiation_on())):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing EM_radiation_on."
	pause(__name__)

	##The two pause functions and show_graph do not need test code.
	
	##Test all_are_numbers:##
	print "\n--------------------------------------------------\n"
	print "Testing all_are_numbers:\n"
	print "Given all numbers returns: " , all_are_numbers(tAllNumbers)
	print "Given not all numbers returns: " , all_are_numbers(tNotAllNumbers)
	tResults = []
	for tNumber in tNotAllNumbers:
		print type(tNumber) , tNumber , "gives: " , all_are_numbers([tNumber])
		tResults.append(all_are_numbers([tNumber]))
	if sum(tResults) == 8:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing all_are_numbers."
	pause(__name__)

	##Test force_float_number:##
	print "\n--------------------------------------------------\n"
	print "Testing force_float_number:\n"
	print "Given an integer returns: " , force_float_number(tAnInteger), type(force_float_number(tAnInteger))
	print "Given a float returns: " , force_float_number(tAFloat), type(force_float_number(tAFloat))
	print "Given a long returns: " , force_float_number(tALong), type(force_float_number(tALong))
	if (type(force_float_number(tAnInteger)) == float) and (type(force_float_number(tAFloat)) == float) and (type(force_float_number(tALong)) == float):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing force_float_number."
	pause(__name__)

	##Test check_float:##
	print "\n--------------------------------------------------\n"
	print "Testing check_float:\n"
	print "Given an integer returns: " , check_float(tAnInteger)
	print "Given a float returns: " , check_float(tAFloat)
	print "Given a long returns: " , check_float(tALong)
	if (not check_float(tAnInteger)) and (check_float(tAFloat)) and (not check_float(tALong)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_float."
	pause(__name__)

	##Test valid_four_vector_index:##
	print "\n--------------------------------------------------\n"
	print "Testing valid_four_vector_index:\n"
	tResults = []
	for tIndex in tIndices:
		print tIndex , "gives: " , valid_four_vector_index(tIndex)
		tResults.append(valid_four_vector_index(tIndex))
	if sum(tResults) == 4:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing valid_four_vector_index."
	pause(__name__)

	##Test valid_dipole_index:##
	print "\n--------------------------------------------------\n"
	print "Testing valid_dipole_index:\n"
	tResults = []
	for tIndex in tIndices:
		print tIndex , "gives: " , valid_dipole_index(tIndex)
		tResults.append(valid_dipole_index(tIndex))
	if sum(tResults) == 2:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing valid_dipole_index."
	pause(__name__)

	##Test valid_particle_status_code:##
	print "\n--------------------------------------------------\n"
	print "Testing valid_particle_status_code:\n"
	tResults = []
	for tCode in tIndices:
		print tCode , "gives: " , valid_particle_status_code(tCode)
		tResults.append(valid_particle_status_code(tCode))
	if sum(tResults) == 3:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing valid_particle_status_code."
	pause(__name__)

	##Test safe_division:##
	print "\n--------------------------------------------------\n"
	print "Testing safe_division:\n"
	tResults = []
	for tDenominator in tDenominators:
		print tDenominator , "gives: " , safe_division(tDenominator)
		tResults.append(safe_division(tDenominator))
	if sum(tResults) == 9:
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing safe_division."
	pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "////////////////////////////////////"
	print "Finished checking assertions module!"
	print "////////////////////////////////////"