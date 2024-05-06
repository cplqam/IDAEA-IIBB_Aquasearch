Aquasearch vs 1.0.0 (version may 6, 2024)

Aquasearch is a tool for the first proteomics characterization of samples analyzed 
using MALDI-TOF in a rapid way. 

Table of contents
-----------------
1.- DEPENDENCIES
2.- INSTALLATION
3.- INPUT FILE FORMAT
4.- CONTRIBUTE


DEPENDENCIES
----------------
If you work with Aquasearch source you need Python installed (tested with Python
version 3.8.10) as well as the following dependencies:
    Matplotlib
    Pandas 
    wxPython 4.2.1
    scikit-learn


INSTALLATION
----------------
 Tested in windows win10 64 bits

Install Aquasearch FROM GitHub (when available)
    1.- Download IDAEA-IIBB_Aquasearch zip archive from repository
    2.- Unzip in python site-packages or in another folder in the python path
    3.- Run Aquasearch by double-clicking on the Aquasearch.pyw module.
    

INPUT FILE FORMAT
-------------------------------------
Aquasearch accepts input files in TXT format for new queries. This files are directly
exported from the MALDI-TOF instrument and have the m/z values 1 one column with 
their intensities in other column

New entry to databese files in mixture samples require different types of samples:
  1. Mixture analyzed by MALDI-TOF (with the same format as new query files)
  2. Mixture analyzed by LC-HRMS and Thermo Proteome Discoverer and exported in XLSX format:
       a) Results of peptide identifications, only the columns with the following
          headers are needed: "A4", "Sequence", "# PSMs" "# Proteins", "# Protein Groups",
          "Protein Group Accessions", "Modifications", "ΔCn", "q-Value" "PEP", "XCorr",
          "Charge", "MH+", "[Da]", "ΔM [ppm]", "RT [min]", "# Missed Cleavages"
       b) Results of protein identifications, only the columns with the following
          headers are needed: "Accession", "Description", "Score", "Coverage", "# Proteins",
          "# Unique Peptides", "# Peptides", "# PSMs", "# AAs", "MW [kDa]", "calc. pI"
  3. Standard analyzed by MALDI-TOF (with the same format as new query files)
  4. Standard unique/non unique peptides obtained using MS-Digest and PeptideRank and exported
     in XLSX format require the following columns: "Number", "m/z (mi)", "m/z (av)", "Charge",
     "Modifications", "Start AA", "End AA", "Missed Cleavages", "Previous AA", "Sequence",
     "Next AA"


CONTRIBUTE
-----------

These programs are made to be useful. If you use them and don't work
entirely as you expect, let me know about features you think are needed, and
I'll try to include them in future releases.

Carlos Perez Lopez may 2024
