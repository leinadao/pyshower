####~~ PyShower 1.0 ~~####
###Copyright 2015/16, Daniel Osborne, All Rights Reserved###
##Durham Thesis: 'Simulations for Particle Physics: Implementing the Colour Dipole Model with Invariant Transverse Momentum Ordering'.##
##For: MPhys Theoretical Physics.##

"""A module setting the precision to be used throughout dipole showering."""

##Deals with the fact that pythons precision limits the precision of some of the showering methods.

##Import required modules:##
import assertions

print "\n/////////////////////////"
print "Loading precision module:"
print "/////////////////////////\n"

##Functions:##

def precise_to():
	"""A function which returns the precision to which the code is being run."""
	__precisionValue = 1.0
	for __i1 in range(constants.number_decimal_places()):
		__precisionValue /= 10.0
	return __precisionValue

def division_precise_to():
	"""A function which returns the precision to which the code is being run when checking for division by zero."""
	__precisionValue = 1.0
	for __i1 in range(constants.division_number_decimal_places()):
		__precisionValue /= 10.0
	return __precisionValue

def round_python_errors(numberToRound):
	"""A function for founding a float at the number of decimal places used throughout."""
	assert assertions.all_are_numbers([numberToRound])
	assert assertions.check_float(numberToRound)
	__numberDecimalPlaces = constants.number_decimal_places()
	return round(numberToRound,__numberDecimalPlaces)

def check_numbers_equal(number1,number2):
	"""A function to check if two numbers are equal to the number of decimal places used throughout."""
	assert assertions.all_are_numbers([number1,number2])
	if (abs(float(number1 - number2)) > precise_to()):
		return False
	return True

def division_check_numbers_equal(number1,number2):
	"""A function to check if two numbers are equal to the number of decimal places used throughout for division."""
	assert assertions.all_are_numbers([number1,number2])
	if (abs(float(number1 - number2)) > division_precise_to()):
		return False
	return True

##Classes:##

##Import here to prevent cyclic import error:
import constants

##Module test code:##
if __name__ == "__main__":
	##Import modules required for testing:##

	##Begin testing:##
	print "\n----------------------------------------------------------------------"
	print "----------------------------------------------------------------------\n"
	print "/////////////////////////"
	print "Testing precision module:"
	print "/////////////////////////"
	assertions.pause(__name__)

	##Setup here:##
	print "\nGenerating test values..."
	tNDP1 = constants.number_decimal_places()
	tNDP2 = constants.division_number_decimal_places()
	tValue1 = 23.2654859661251853
	tValue2 = 1.0 * tValue1
	tRoundedTemp = "23.2654860"
	tRoundedShouldBe = float(tRoundedTemp[:tNDP1+3])

	##Test precise_to:##
	print "\n--------------------------------------------------\n"
	print "Testing precise_to:\n"
	print "Calling number_decimal_places returns: " , constants.number_decimal_places()
	print "Calling precise_to returns: " , precise_to()
	if (check_numbers_equal(precise_to(), 0.000000001)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print " \nFinished testing precise_to."
	assertions.pause(__name__)

	##Test division_precise_to:##
	print "\n--------------------------------------------------\n"
	print "Testing division_precise_to:\n"
	print "Calling number_decimal_places returns: " , constants.division_number_decimal_places()
	print "Calling division_precise_to returns: " , division_precise_to()
	if (check_numbers_equal(division_precise_to(), 0.00000000001)):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print " \nFinished testing division_precise_to."
	assertions.pause(__name__)

	##Test round_python_errors:##
	print "\n--------------------------------------------------\n"
	print "Testing round_python_errors:\n"
	print "Given the test value: {0:.20f}".format(tValue1)
	tResult = round_python_errors(tValue1)
	print "Returns the rounded result: {0:.20f}".format(tResult)
	print "Printing it without formatting gives:", tResult
	if (tRoundedShouldBe == tResult):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
		print "\nPerhaps tRoundedTemp in testCode hasn't been updated?"
	print "\nFinished testing round_python_errors."
	assertions.pause(__name__)

	##Test check_numbers_equal:##
	print "\n--------------------------------------------------\n"
	print "Testing check_numbers_equal:\n"
	print "Using a test number of:" , tValue1
	print "With:", tValue2 , "returns:" , check_numbers_equal(tValue1,tValue2)
	print "Adding", precise_to()/10, "to the second value gives:" , check_numbers_equal(tValue1,(tValue2 + precise_to()/10))
	print "Adding", precise_to(), "to the second value gives:" , check_numbers_equal(tValue1,(tValue2 + precise_to()))
	print "Adding", precise_to()*10, "to the second value gives:" , check_numbers_equal(tValue1,(tValue2 + precise_to()*10))
	tResults = [check_numbers_equal(tValue1,(tValue2 + precise_to()/10)),check_numbers_equal(tValue1,(tValue2 + precise_to()))]
	tResults += [check_numbers_equal(tValue1,(tValue2 + precise_to()*10))]
	if (sum(tResults) == 1):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing check_numbers_equal."
	assertions.pause(__name__)

	##Test division_check_numbers_equal:##
	print "\n--------------------------------------------------\n"
	print "Testing division_check_numbers_equal:\n"
	print "Using a test number of:" , tValue1
	print "With:", tValue2 , "returns:" , division_check_numbers_equal(tValue1,tValue2)
	print "Adding", division_precise_to()/10, "to the second value gives:" , division_check_numbers_equal(tValue1,(tValue2 + division_precise_to()/10))
	print "Adding", division_precise_to(), "to the second value gives:" , division_check_numbers_equal(tValue1,(tValue2 + division_precise_to()))
	print "Adding", division_precise_to()*10, "to the second value gives:" , division_check_numbers_equal(tValue1,(tValue2 + division_precise_to()*10))
	tResults = [division_check_numbers_equal(tValue1,(tValue2 + division_precise_to()/10))]
	tResults += [division_check_numbers_equal(tValue1,(tValue2 + division_precise_to()))]
	tResults += [division_check_numbers_equal(tValue1,(tValue2 + division_precise_to()*10))]
	if (sum(tResults) == 1):
		print "\nTest successful!"
	else:
		print "\nTest failed!"
	print "\nFinished testing division_check_numbers_equal."
	assertions.pause(__name__)

	##Done testing:##
	print "\n---------------------------------------------\n"
	print "///////////////////////////////////"
	print "Finished checking precision module!"
	print "///////////////////////////////////"