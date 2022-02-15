P1812 Version 2 (28.05.18)

MATLAB implementation of Recommendation ITU-R P.1812-4

GENERAL NOTES
--------------

Files and subfolders in the distribution .zip package.

 tl_p1812.m                  - MATLAB function implementing Recommendation ITU-R P.1812-4.

 validate_p1812.m            - MATLAB script used to validate the implementation of 
                               Recommendation ITU-R P.1812-4 as defined in the file 
                               tl_p1812.m using a set of test terrain profiles provided 
                               in the folder ./validation_profiles/

 ./validation_profiles/	     - Folder containing a proposed set of terrain profiles for
                               validation of MATLAB implementation (or any other software
                               implementation) of Recommendation ITU-R P.1812-4

 validateion_result_log.csv  - Template for reporting final and intermediate results of basic
                               transmission loss computation according to Recommendation ITU-R P.1546-5

 ./validation_results/       - Folder containing all the results written during the transmission loss
                               computations for the set of terrain profiles defined
                               in the folder ./validation_profiles/

 ./src/                      - Folder containing the functions used by tl_p1812.m and validate_p1812.m


UPDATES AND FIXES
-----------------
Version 2 (28.05.18)

	- Corrected a bug in printing the values for Lbulla50 and Lbulls50 (scalar and not vector values) reported out by Damian Bevan
    - Corrected a typo in error reporting for hrx
	- Introduced path center calculation according to Annex H of Recommendation ITU-R P.2001-2
	- Re-created validation examples so that they match the new path-center calculation (minor differences of the order 0.01 dB)
	- Created two new validation examples to test for the case of vertical polarization
    - Added function isOctave to ensure compatibility with interactive plotting in Octave



Version 1 (18.11.16)
     - Initial implementation
