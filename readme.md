Binary Executables
------------------
- lbc_smokehist_win64.exe: compiled on Cygwin on Windows 7 64-bit using GCC version 3.4.4 and Make 3.81.
- lbc_smokehist_ubuntu64.exe: compiled on Ubuntu 12.04 64-bit using GCC compiler version 4.6.3 and Make 3.81.
- lbc_smokehist_osx64.exe: compiled on Mac OS X 64-bit version 10.7.5 using GCC 4.2.1 and Make 3.81.

Compiling from Source
---------------------
Windows systems:
- The use of Cygwin is the recommended means of compiling the second generation shg on windows systems.
- Cygwin may be obtained by downloading and installing the executable: http://cygwin.com/setup.exe.
- By accepting all defaults and selecting a local mirror location, the default installation package appropriately installs GCC v3.4.4 and Make v3.8.1 as listed in the compatability section above.

All systems:
- Issue the command `sh install.sh` from the project root (more details found in the makefile).
- Note: It may be necessary to ensure the install.sh file is executable with the command `chmod +x install.sh`.
- An executable file named lbc_smokehist.exe (by default) should be created in the project root.

Quick Start
-----------
- Following compilation, a version 6.2.0 test run with a cutoff year of 2050 could be performed from the project root using:
`./lbc_smokehist.exe data/shg2p0 1 2 3 4 test.in test.out 1 0 -c 2050`
- Note, if using a binary, the lbc_smokehist_*.exe selected needs to correspond to your OS.
- See the HelpFile.txt for more details regarding details of command line usage.
- Running the python tests may be accomplished with `python run_tests.py 1 > test.out`. This has a numpy dependency. I defer to a google search for the proper setup of numpy for python on your OS.

Questions
---------
Please email any concerns to ben.racine@cornerstonenw.com
Any modifications will result in a minor release and be updated accordingly on the SHG website.
