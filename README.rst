############################
T4ME - Transport 4 MatErials
############################

       _________________  ____  ___    _______________

      /            /   / /   / /   \__/   /          /

     /____    ____/   / /   / /          /   _______/ 

         /   /   /   /_/   /_/          /   /___      

        /   /   /           /   /\_/   /   ____/      

       /   /   /_____    __/   /  /   /   /_______    

      /   /         /   / /   /  /   /           /    

     /___/         /___/ /___/  /___/___________/    


Routines to calculate the transport properties of materials
using the linearized Boltzmann Transport Equations (BTE)
in the Relaxtion-Time-Approximation (RTA).

Features
********

Here comes a condensed list of features.

Structure
*********

The structure of the program is simple: the main routines
are written in Python utlizing NumPy and SciPy where
necessary. In addition there are calls to external
routines through Cython, particularly the optional libraries.
It should support Python > 2.7 and Python 3.

The main driver for the program is the :file:`t4me.py` file
located in the main directory. This is used to execute
the program and call necessary subroutines.

Contributing and versioning
***************************

Standard Git versioning is utilized. Contributions are welcome,
encouraged and (greatly) appreciated.

Author
******

Espen Flage-Larsen with finances from the Norwegian
Research Council, Thelma project (228854).

License
*******

This project is licensed under the GNU GPLv3. Please see
:file:`LICENSE.md` for additional details.
