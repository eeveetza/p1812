P1812 Version 6.0 (22.05.22)

MATLAB implementation of Recommendation ITU-R P.1812-6

GENERAL NOTES
--------------

Files and subfolders in the distribution .zip package.

 tl_p1812.m                  - MATLAB function implementing Recommendation ITU-R P.1812-6.

 validate_p1812.m            - MATLAB script used to validate the implementation of 
                               Recommendation ITU-R P.1812-6 as defined in the file 
                               tl_p1812.m using a set of test terrain profiles provided 
                               in the folder ./validation_profiles/

 ./validation_profiles/	     - Folder containing a proposed set of terrain profiles for
                               validation of MATLAB implementation (or any other software
                               implementation) of Recommendation ITU-R P.1812-6

 validation_result_log.csv  - Template for reporting final and intermediate results of basic
                               transmission loss computation according to Recommendation ITU-R P.1812-6

 ./validation_results/       - Folder containing all the results written during the transmission loss
                               computations for the set of terrain profiles defined
                               in the folder ./validation_profiles/

 ./private/                  - Folder containing the functions used by tl_p1812.m and validate_p1812.m


UPDATES AND FIXES
-----------------

Version 6.0 (22.05.22)

	- Aligned to ITU-R P.1812-6 (excluded terminal clutter model, aligned free-space loss, extended frequency limit)
	- Renamed subfolder "src" into "private" which is automatically in the MATLAB search path
	- Simplified handling of (optional) input arguments using name-value pairs 
	- Ensured that the variable series is a row vector in find_intervals.m
	- Used 2.998e8 m/s for speed of light as per ITU-R P.2001-4 (instead of 3e8 m/s)
        - Use comma instead of semicolon as a separator in the output csv files 

Version 5.4 (08.04.21)

	- Introduced location variability for outdoor propagation
        - Introduced alternative method to compute Lbulls w/o using 
          terrain profile (Attachment 4 to Annex 1)

Version 4.1 (28.05.18)

	- Corrected a bug (swap of long and lat coordinates in great_circle_path call, reported by Wladislaw Budynkiewicz.
	- Corrected a bug in printing the values for Lbulla50 and Lbulls50 (scalar and not vector values) reported out by Damian Bevan
        - Corrected a typo in error reporting for hrx
	- Introduced path center calculation according to Annex H of Recommendation ITU-R P.2001-2
	- Re-created validation examples so that they match the new path-center calculation (minor differences of the order 0.01 dB)
	- Created two new validation examples to test for the case of vertical polarization
        - Added function isOctave to ensure compatibility with interactive plotting in Octave
        - Added an error message occuring in the case ./src/ is not on the MATLAB search path



Version 1 (18.11.16)
     - Initial implementation


License and copyright notice

Swiss Federal Office of Communications OFCOM (hereinafter the "Software Copyright Holder") makes the accompanying software 
(hereinafter the "Software") available free from copyright restriction. 

The Software Copyright Holder represents and warrants that to the best of its knowledge, 
it has the necessary copyright rights to waive all of the copyright rights as permissible under national law in the Software 
such that the Software can be used by implementers without further licensing concerns. 

No patent licence is granted, nor is a patent licensing commitment made, by implication, estoppel or otherwise. 

Disclaimer: Other than as expressly provided herein, 

(1) the Software is provided "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO, 
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGMENT OF INTELLECTUAL PROPERTY RIGHTS and 

(2) neither the Software Copyright Holder (or its affiliates) nor the ITU shall be held liable in any event for any damages whatsoever 
(including, without limitation, damages for loss of profits, business interruption, loss of information, or any other pecuniary loss) 
arising out of or related to the use of or inability to use the Software.