EXE = matrix_lib_test
ADA_VERSION = -gnat05
SRC = Number_Theory_Tools.adb Permutations_Generic.adb Matrix_Lib.adb Matrix_Lib_Test.adb 
FLAGS = -gnato -gnatwa -fstack-check

all: 
	gnatmake $(ADA_VERSION) $(FLAGS) $(SRC)

ada83: 
	gnatmake -gnat83  $(FLAGS) $(SRC)

ada95: 
	gnatmake -gnat95  $(FLAGS) $(SRC)

ada2005: 
	gnatmake -gnat05  $(FLAGS) $(SRC)

ada2012: 
	gnatmake -gnat12  $(FLAGS) $(SRC)

clean:
	rm *.ali *~ *.o $(EXE)