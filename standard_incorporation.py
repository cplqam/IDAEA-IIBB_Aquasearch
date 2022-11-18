import pandas as pd
import numpy
import os
import protein_signals as ps
import sqlite_requests as sr

def xlsx_extraction(path_):
    """ This function relates each peptide and its unique with its mz value
    
    INPUT
    path_: string. It is the path of the query file
    
    
    >>> result = xlsx_extraction('test_files/Standares/Aquasearch_Proteins_Unique.xlsx')
    >>> result.shape[1] == 5
    True
    
    """
    
    df_unique = pd.read_excel(path_, sheet_name=0) 
    protein_codes = df_unique['protein'].tolist()
    protein_codes_filtered = []
    
    for code in protein_codes:
        code = code.split('|')
        protein_codes_filtered.append(code[1])
        
    df_unique['protein'] = protein_codes_filtered  
    df_unique = df_unique.loc[:,['protein', 'peptide', 'unique']]
    
    
    df_mz = pd.read_excel(path_, sheet_name = 2) 
    df_mz = df_mz.loc[:,['Number', 'm/z (mi)', 'm/z (av)', 'Sequence']]

    list_code = []
    list_mz = []
    list_pep = []
    list_uniq = []
    list_mm = []

    for code, pep, uniq in zip(df_unique['protein'], df_unique['peptide'], df_unique['unique']):
        for code2, mz, mm, pep2 in zip(df_mz['Number'], df_mz['m/z (mi)'], df_mz['m/z (av)'], df_mz['Sequence']):
            
            if (code == code2) & (pep == pep2):
                list_code.append(code)
                list_pep.append(pep)
                list_mz.append(mz)
                list_mm.append(mm)
                if uniq == 'unique':
                    list_uniq.append(1)
                elif uniq == 'ambiguous':
                    list_uniq.append(0)
    df_final = pd.DataFrame({'Code': list_code,'Peptide': list_pep, 'mz': list_mz, 'Unique': list_uniq, 'molecular_mass': list_mm}) 
    return df_final

def import_txt(path_):
    """This function import a table from txt file
        
        >>> df_2, mzint_2 = import_txt('test_files/Standares/pmf_H1.txt')
        >>> df_2.shape == mzint_2.shape
        True
    
    """
    
    df = pd.read_csv(path_, sep = '\s+')
    df = df.iloc[:,[0,1]]
    
    df.columns =  ['mz', 'intensity']
    
    mz_s = df['mz'].squeeze()
    int_s = df['intensity'].squeeze()

    mz_int = numpy.array([mz_s, int_s]).transpose()
    
    return(df, mz_int)



def standard_complete(maldi, path_, code, sample, ppm = 100, db = 'Aquasearch_study'):
    """This function completes the unique/non unique information for each signal obtained from
       standard analysis by MALD-TOF. This function works with the auxiliar function 
       xlsx_extraction()
       
       INPUT
       maldi: string. It is the path of the maldi file which needs to be completed
       path_: string. It is the path of the file whit the unique/non unique info
       code: string. The uniprot code of the protein standad analyzed by MALDI-TOF
       ppm: integer. Maximum error allowed for considering signals as the same.
            By default --> ppm = 100
    """

    df_unique = xlsx_extraction(path_)
    df_unique = df_unique[df_unique['Code'] == code]
    
    df_maldi, array_maldi = import_txt(maldi)
    
    uniq_inf = []
    seq_inf = []
    mm_inf = []
    stand_inf = []
    sampl_inf = []
    
    for i in range(df_maldi.shape[0]):
        mz_query = df_maldi.loc[i, 'mz']
        coincidence = 0
    
        for j in range(df_unique.shape[0]):
            mz_ans = df_unique.iloc[j, 2]
            uniq = df_unique.iloc[j, 3]
            seq = df_unique.iloc[j, 1]
            mm = df_unique.iloc[j, 4]
            if mz_ans > mz_query:
                ppm_calculated = (1 - (mz_query / mz_ans)) * 1E6
            else:
                ppm_calculated = (1 - (mz_ans / mz_query)) * 1E6

            if ppm_calculated <= ppm:
                coincidence += 1
                uniq_val = uniq
                seq_val = seq
                mm_val = mm
        
        if coincidence > 0:
            uniq_inf.append(str(uniq_val))
            stand_inf.append('1')
            seq_inf.append(seq_val)
            mm_inf.append(mm_val)
            sampl_inf.append(sample)
        else:
            uniq_inf.append('-')
            stand_inf.append('-')
            seq_inf.append('-')
            mm_inf.append('-')
            sampl_inf.append('-')
                
    df_maldi['Unique Pep'] = uniq_inf
    df_maldi['Standard signal'] = stand_inf
    df_maldi['Peptide seq'] = seq_inf
    df_maldi['Peptide mass'] = mm_inf
    df_maldi['Sample'] = sampl_inf
        
    df_maldi = df_maldi[(df_maldi['Unique Pep'] != '-')]
    df_maldi = ps.relat_intensity_calc(df_maldi)
    
    exp = df_maldi[['mz', 'intensity', 'Unique Pep']]
    maldi = maldi.replace('.txt','') 
    maldi_completed = maldi + '_completed.txt'
    exp.to_csv(maldi_completed, header='Project', index=None, sep= '\t')
  
    id_p = ps.table_proteins_check(code)
    code_seq = []
    for seq, uniq, mass in zip(df_maldi.loc[:, 'Peptide seq'], df_maldi.loc[:, 'Unique Pep'], df_maldi.loc[:, 'Peptide mass']):
        id_s = ps.table_sequence_check(seq, uniq, mass)
        code_seq.append(id_s)
    
    df_maldi.loc[:, 'Peptide seq'] = code_seq
    
    table_examined = ps.table_spectrum_check(df_maldi, id_p)
    sr.eliminate_table(db, 'Spectrums_table')
    sr.spectrums_table(db)
    sr.insert_new_spectrum(db, table_examined)
    
    ps.quality_filter(db)
    ps.reference_spectrum(db)
    
    
def standard_complete_union(maldi, path_, code_prot, code_stand, sample, ppm = 100, db = 'Aquasearch_study'):
    """This function completes the unique/non unique information for each signal obtained from
       standard analysis by MALD-TOF. This function works with the auxiliar function 
       xlsx_extraction()
       
       INPUT
       maldi: string. It is the path of the maldi file which needs to be completed
       path_: string. It is the path of the file whit the unique/non unique info
       code_prot: string. The uniprot code of the protein standad analyzed by MALDI-TOF
       code_stand: string. The name of the new estry resulting fron the union of 2 proteins
       sample: string. The name of the sample
       ppm: integer. Maximum error allowed for considering signals as the same.
            By default --> ppm = 100
       db: string. The name of the database
    """
    df_unique = xlsx_extraction(path_)
    df_unique = df_unique[df_unique['Code'] == code_prot]

    df_maldi, array_maldi = import_txt(maldi)


    uniq_inf = []
    seq_inf = []
    mm_inf = []
    stand_inf = []
    sampl_inf = []

    for i in range(df_maldi.shape[0]):
        mz_query = df_maldi.loc[i, 'mz']
        coincidence = 0

        for j in range(df_unique.shape[0]):
            mz_ans = df_unique.iloc[j, 2]
            uniq = df_unique.iloc[j, 3]
            seq = df_unique.iloc[j, 1]
            mm = df_unique.iloc[j, 4]
            if mz_ans > mz_query:
                ppm_calculated = (1 - (mz_query / mz_ans)) * 1E6
            else:
                ppm_calculated = (1 - (mz_ans / mz_query)) * 1E6

            if ppm_calculated <= ppm:
                coincidence += 1
                uniq_val = uniq
                seq_val = seq
                mm_val = mm
        
        if coincidence > 0:
            uniq_inf.append(str(uniq_val))
            stand_inf.append('1')
            seq_inf.append(seq_val)
            mm_inf.append(mm_val)
            sampl_inf.append(sample)
        else:
            uniq_inf.append('-')
            stand_inf.append('-')
            seq_inf.append('-')
            mm_inf.append('-')
            sampl_inf.append('-')
                
    df_maldi['Unique Pep'] = uniq_inf
    df_maldi['Standard signal'] = stand_inf
    df_maldi['Peptide seq'] = seq_inf
    df_maldi['Peptide mass'] = mm_inf
    df_maldi['Sample'] = sampl_inf
        
    df_maldi = df_maldi[(df_maldi['Unique Pep'] != '-')]
    df_maldi = ps.relat_intensity_calc(df_maldi)

    exp = df_maldi[['mz', 'intensity', 'Unique Pep']]
    maldi = maldi.replace('.txt','') 
    maldi_completed = maldi + '_completed.txt'
    exp.to_csv(maldi_completed, header='Project', index=None, sep= '\t')

    id_p = ps.table_proteins_check(code_stand)
    code_seq = []
    for seq, uniq, mass in zip(df_maldi.loc[:, 'Peptide seq'], df_maldi.loc[:, 'Unique Pep'], df_maldi.loc[:, 'Peptide mass']):
        id_s = ps.table_sequence_check(seq, uniq, mass)
        code_seq.append(id_s)

    df_maldi.loc[:, 'Peptide seq'] = code_seq

    table_examined = ps.table_spectrum_check(df_maldi, id_p)
    sr.eliminate_table(db, 'Spectrums_table')
    sr.spectrums_table(db)
    sr.insert_new_spectrum(db, table_examined)

    ps.quality_filter(db)
    ps.reference_spectrum(db)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
    
    
    path__ = 'test_files/Standares'
    
    unique_file = 'Aquasearch_Proteins_Unique.xlsx'
    maldi_file = 'pmf_H10.txt'
    
    unique_path = os.path.join(path__, unique_file)
    maldi_path = os.path.join(path__, maldi_file)
    
    df_unique = xlsx_extraction(unique_path)
    
    standard_complete(maldi_path, unique_path,'P08835', 'standard_1', ppm = 100)
    