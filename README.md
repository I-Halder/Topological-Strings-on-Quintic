# **Topological-Strings-on-Quintic**

## Folder P
This folder contains information on the free energy of topological strings on quintic up to genus 49. To use it, follow the steps below
  1. Download the files in a local folder called "P". In MacOS, open the terminal and access the folder "P"
  2. Use the following command to obtain a Mathematica readable file named P.mx: **cat * > P.mx**

## File ABGt.mx
This file contains information on the moduli space geometry and Yamaguchi-Yau variables.

## File example.nb
Put the files example.nb, P.mx, ABGt.mx in the same folder. Now open example.nb in Mathematica and follow the instructions.

## File gv.csv
This file contains Gopakumar-Vafa coefficients for quintic. Each column represents data for a given genus (0 to 49), and each row corresponds to a particular degree (1 to 21).

## File quintic.m
This is the actual Mathematica file we used for the calculation. An initial version of it was written by IH, and the current optimized form is due to YL.

Customize genus, degree, GUI/command line, and directory in Initialization/Input cell. For command line, run e.g. "math -script quintic.m -g 7" where math is alias for Mathematica executable.
