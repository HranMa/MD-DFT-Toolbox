# md-scripts

This project is designed to store some auxiliary scripts written with ChatGPT during MD simulation.
(Of course, DFT calculations including, too)

Most of the scripts can be run in the WSL with Ubuntu 22.04 installed, since it's much convenient to set up than normal Windows enviroment. 

If you're a graduate student under the guidence of Prof. Wang Dan or you have used the scripts in this project, remember to Star it.


If you have any further questions, please contact me through this email: mahaorancn@163.com

## ice-growth-front
A series scripts to determine the ice growth front for Gromacs simulation results. 
To be more precisely, these scripts can be used in the simulation boxes with the ice in the middle and water on the upper and lower two sides. 
For ResNames, ice == fSOL and water == SOL.
TIP4P/ice was used, thus in the .gro files, water molecules have 4 atoms.

## msd
A series scripts to determine some ions positions in the simualtion. 

## normal-tools
A toolkit caters to daily needs
