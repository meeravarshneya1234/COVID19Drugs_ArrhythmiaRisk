![Title](https://github.com/meeravarshneya1234/COVID19Drugs_ArrhythmiaRisk/blob/master/Images/Title.png)  

This folder contains source code to reproduce the figures for "Investigational treatments for COVID-19 may increase ventricular arrhythmia risk through drug interactions" by Varshneya & Irurzun-Arana et. al. This work was performed in the Cardiac Systems Pharmacology Lab of Dr. Eric A. Sobie in the Department of Pharmacological Sciences at the Icahn School of Medicine Mount Sinai. 

## Requirements
MATLAB - version 2014 or higher; 2017a was used to create the populations and 2019b was used to perform analysis.

R - version 3.6.3; Required packages: mrgsolve, ggplot2, dplyr, Rtools

## Usage 
Figures can be recreated by running the "FigureX.m" files. Functions used within each m-file are within the functions folder and have been extensively commented for easy usability.

## Models 
| PK Model        | QSP Model         | 
| :---:         |     :---:      |
| [**Lopinavir + Ritonavir Model**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3962720/) | [**O'Hara Ventricular Myocyte Model**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002061) | 
| [**Azithromycin Model**](https://aac.asm.org/content/58/11/6675.long) |
| [**Chloroquine Model**](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/cpt.1893) |

## Abstract
Many drugs that have been proposed for treatment of COVID-19 are reported to cause cardiac adverse events, including ventricular arrhythmias. In order to properly weigh risks against potential benefits, particularly when decisions must be made quickly, mathematical modeling of both drug disposition and drug action can be useful for predicting patient response and making informed decisions. Here we explored the potential effects on cardiac electrophysiology of 4 drugs proposed to treat COVID-19: lopinavir, ritonavir, chloroquine, and azithromycin, as well as combination therapy involving these drugs. Our study combined simulations of pharmacokinetics (PK) with quantitative systems pharmacology (QSP) modeling of ventricular myocytes to predict potential cardiac adverse events caused by these treatments. Simulation results predicted that drug combinations can lead to greater cellular action potential prolongation, analogous to QT prolongation, compared with drugs given in isolation. The combination effect can result from both pharmacokinetic and pharmacodynamic drug interactions. Importantly, simulations of different patient groups predicted that females with pre-existing heart disease are especially susceptible to drug-induced arrhythmias, compared with diseased males or healthy individuals of either sex. Statistical analysis of population simulations revealed the molecular factors that certain females with heart failure especially susceptible to arrhythmias. Overall, the results illustrate how PK and QSP modeling may be combined to more precisely predict cardiac risks of COVID-19 therapies. 

## Questions/Contributions
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Citation 
Varshneya, M., Irurzun-Arana, I., Campana, C., Dariolli, R., Gutierrez, A., Pullinger, T.K. and Sobie, E.A. (2021), Investigational Treatments for COVID‐19 may Increase Ventricular Arrhythmia Risk Through Drug Interactions. CPT Pharmacometrics Syst. Pharmacol., 10: 100-107. https://doi.org/10.1002/psp4.12573
