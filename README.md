Vicsek 2D Simulation
====================

Steps to compile
----------------

At the very begining of each .f90 file in the src folder you will find a line like this:

> !COMPILAMOS
> ! \> g95 ....

copy the second line (skiping the ! \> symbol) into your terminal.

1. Compile the modules 'ModuloSimulaciones.f90' and 'ModuloSim_Vicsek.f90' to get the objects for the executable.

2. Compile the executables 'Vicsek.f90' and 'analisisACTc.f90'

Steps to run
------------

1. Fill the InputFile 'Entrada_Parametros_Vicsek' with the desired initial parameters.

2. Execute 'Vicsek.exe': choosing the noise-model and other choises that will pop up.
##3rd 

For the analysis check the Vicsek_Implementation file
-----------------------------------------------------