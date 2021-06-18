# Vicsek 2D Simulation

![](docs/PosterRAFA_Simulaciones-page-001.jpg)
## Overview

This repo contains both the code for the Vicsek2D simulation and the presentation docs (in spanish),
which includes the analysis report and the poster.\
This project aimed to understand the impact of two different dynamic rules over the simulation of an autonomous system,
in the context of Statistical Mechanics under the Vicsek model hypothesis.
In particular,
we explore the impact of the dynamic with regards to the order parameter to determine whether the phase transition is the same for both dynamics.
We conclude that they indeed undergo a 1st order transition,
and we estimate the critical points for each one.


## Steps to compile

At the very begining of each .f90 file in the src folder you will find a line like this:

> !COMPILAMOS

> ! \> g95 ....

copy the second line (skiping the ! \> symbol) into your terminal.

1. Compile the modules '**ModuloSimulaciones.f90**' and '**ModuloSim_Vicsek.f90**' to get the objects for the executable.

2. Compile the executables '**Vicsek.f90**' and '**analisisACTc.f90**'

## Steps to run

1. Fill the InputFile '**Entrada_Parametros_Vicsek**' with the desired initial parameters.

2. Run '**Vicsek2D.exe**': choosing the noise-model and other choises that will pop up.

3. Run '**AnalisisACTc**'

## For the analysis check the Vicsek_Implementation file

Full notes on the logic of the scripts can be found in the **Memorias Vicsek** file.
