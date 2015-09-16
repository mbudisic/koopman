# KOOPMAN MODE DECOMPOSITION #
Matlab toolbox for Koopman mode decomposition.

Written by: Marko Budisic

The goal of this toolbox is to collect several common Koopman mode decomposition algorithms, in a documented, transparent code.

The code is licensed under a BSD3 license, found in the accompanying LICENSE file.

For the content of the package, please consult +koopman/Contents.m.

# Installation

Place the toolbox on the drive, e.g., to ~/MatlabToolbox/koopman, and add its top folder to Matlab path.

    >> addpath("~/MatlabToolbox/koopman")
    >> savepath

Then functions from the toolbox can be accessed via namespace "koopman", e.g.,

    >> koopman.DMD


# Use

Currently, the toolbox implements two algorithms based on Dynamic Mode Decomposition, and one algorithm based on Discrete Fourier Transform. All functions have similar syntax. For further documentation, see help lines for individual functions in the koopman namespace.

A demo (and validation) for the toolbox is located in koopman/validate folder. Once the namespace +koopman is in the path, demo can be run by

    >> cd [toolboxfolder]/validate
    >> DemoKoopmanModes

Demo implements an exponentially growing mode used by Duke et al. Spatial shape of the data set is fixed in demo, while the time behavior is set via arguments to the demo function.

Here is an example of the demo run for decay rate -0.1 and angular frequency 21:

![Data for -0.1+21i set](img/data-01_21i.png "Visualization of the input for time frequency -0.1+i21")
![Results for -0.1+21i set](img/results-01_21i.png ""Visualization of the output for time frequency -0.1+i21"")

Here is an example of the demo run for decay rate 0 and angular frequency 20:

![Data for 0+20i set](img/data-0_20i.png "Visualization of the input for time frequency i20")
![Results for 0+20i set](img/results-0_20i.png "Visualization of the output for time frequency i20")