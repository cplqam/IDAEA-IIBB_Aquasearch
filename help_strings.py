tolerance = """The m/z error in ppm to match the signals from new query and the signals in database

By default: ppm=100"""

database = 'Choose the databese in which the signals of the new query have to be searched'

advanced_details = """These parameters affect the importance of the m/z signal and its intensities
to canculate the scores. If these parameters are modilied, you have to
recalculate the score thresholds to confirm the presence of a protein 
in the sample. 
                                 
By default: 
- m/z power = 1
- Intensity power = 0 
"""

database_ne = 'Choose the databese you want to complete with new signals or proteins'

prot_results = """A list with the name of the analized samples is displayed. Clicking on them 
a new window with the results of that sample is displayed. This new 
window shows the proteins in the chosen database, its name, organism, 
the score of the proteins, the number of used peptides and how many 
of them are uniques
Scores:
    - <4: The presence of the protein is unlikely in the sample
    - 4-5: The presence of the protein is likely in the sample
    - >5: The presence of the protein is sure in the sample

The "protein code" column can be clicked to open a window with the 
peptides found in the sample for that protein and their ppm error    
"""

pcs = 'Select the number of Principal Components in PCA model'

pcx = 'The PC number to show in axis X'

pcy = 'The PC number to show in axis Y'

pcz = 'The PC number to show in axis Z'

title = 'Write the title of the PCA figure'

size = 'Choose the size of the points in the PCA figure'

maldi_ind_mix = """Select a MALDI-TOF spectrum in txt format to introduce new signals in the
 selected database"""

pep_ind_mix = """Select a Proteome Discoverer peptide results file in xlsx format 
to introduce new signals in the selected database"""

prot_ind_mix = """Select a Proteome Discoverer protein results file in xlsx format 
to introduce new signals in the selected database"""

ppm_ind_mix = """Select the error (in ppm) to consider a MALDI-TOF signal and an 
identification from Proteome Discoverer as the same signal

By default: ppm=100"""

rank_ind_mix = """Select the number of proteins in abundance of the Proteome Discoverer 
protein results file are considered

By default: Protein rank = 200"""

unique_ind_mix = """Select the identification of the signals related with more than 1 peptide: 
    - Selection: 1 peptide is chosen among all posibilities
    - All posibilities: all the posible peptides are considered
    
By default in Aquasearch_study database: All posibilities"""
    
code_ind_mix = """Type the uniprot code to be searched in the new MALDI-TOF sample 

Ej: P19121
If the option 'all individual protein codes' is selected, all the
individual proteins in database will be searched in the sample and 
this section could be blank or will be ignored"""

name_ind_mix = "The name of the new MALDI-TOF sample added to database"

maldi_indiv_std = """Select a MALDI-TOF standard spectrum in txt format to introduce new 
signals in the selected database"""

unique_ind_std = """Select a xlsx file with the information of the unique/no unique peptide 
in the forst sheet and the information of the peptide fragmentation en the second sheet.

- The information of unique/no unique can be obtained from PeptideRank, 
using as input the fasta files of all proteins in database
- The information of the peptide fragmentation can be obtained from 
Protein Prospector (MS-Digest) using the uniprot protein codes
"""

ppm_ind_std = """Select the error (in ppm) to consider a MALDI-TOF signal and the 
fragmentation from Protein Prospector as the same signal

By default: ppm=100"""

code_ind_std = "Type the uniprot code of the standard"

sample_ind_std = "The name of the new MALDI-TOF standard sample added to database"

code1_gru_mix = """Type the uniprot codes (uniprot code 1 and 2) to be searched and grouped 
in the new MALDI-TOF sample 

Ej: code 1 = P07724 and code 2 = P02770"""

name_gru_mix = """Type the name of the grouped proteins

Ej: Murid albumin
If the option 'all grouped protein codes' is selected, all the grouped proteins
in database will be searched in the sample and this section could be blank 
or will be ignored
"""

code_gru_std = """Type the uniprot code to be searched in the new MALDI-TOF sample

Ej: P19121
"""

name_gru_std = """Type the name of grouped proteins in which the standard peptide 
signals have to be included

Ej: Murid albumin
"""

del_db1 = "Select the database from the protein will be removed"

del_db2 = "Select the database from the sample will be removed"

del_prot = "Select protein to remove"

del_sam = "Select sample to remove"


