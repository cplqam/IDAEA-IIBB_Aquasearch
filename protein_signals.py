import sqlite_requests as sr
import sqlite3
import os
import pandas as pd


def table_comprobation(protein_code, db='Aquasearch_study'):
    """This function checks if the database and a table with the 
       protein accession code name exist. If any of them do not exist,
       it is created
       
       protein_code: string. The name of the table (protein accession code)
       db: string. The name of the DB. By default: Aquasearch_study
       """
    
    # Try to create a database if it doesn't exist
    if not os.path.exists(db):
        sr.create_db(db)
    
    # Try to create a table if it doesn't exist
    try:
        sr.create_table_protein(db, protein_code)
    except sqlite3.OperationalError:
        pass
    

def fill_table(protein_code, maldi_complete, db='Aquasearch_study', options = 1, choose = 1):
    """This function creates or complete a table with the peptide signals belonging to
       the same protein in a mix sample 
       
       protein_code: string. The name of the table (protein accession code)
       maldi_complete: Dataframe. The signals from MALDI with all peptide information.
       Result from pd_maldi_match.xml_complete()
       db: string. The name of the DB. By default: Aquasearch_study
       
       options: integer
       options = 0: there is only 1 peptide option for each signal (pd_maldi_match.xml_complete(unique_ = 0/1))
       options = 1: each signal can have 1 or more peptide signals associated (pd_maldi_match.xml_complete(unique_ = 2))
       
       choose: integer.
       choose = 0: if a new signal is considered as the same as an database signal (according to ppm parameter)
                   only the signal with more relative intensity are considered 
       choose = 1: all the signals are introduced in the database  
    """
    protein_codes = maldi_complete['Protein Accession code']
    protein_codes = protein_codes.tolist()

    if options == 0:
        idx = []
        for i in range(len(protein_codes)):
            if protein_codes[i] == protein_code:
                idx.append(i)
        
        signals_inter = maldi_complete.iloc[idx, :]
        signals_inter = signals_inter.reset_index().drop(['index'], axis=1)
        if signals_inter.shape[0] > 0:
            signals_inter = relat_intensity_calc(signals_inter)
        
    elif options == 1:
        idx = []
        for i in range(len(protein_codes)):
            prot_query = protein_codes[i]
            prot_query = prot_query.split(';')
            
            for prot in prot_query:
                if prot == protein_code:
                    idx.append(i)
                    
        signals_inter = maldi_complete.iloc[idx, :]
        signals_inter = signals_inter.reset_index().drop(['index'], axis=1)
        if signals_inter.shape[0] > 0:
            signals_inter = relat_intensity_calc(signals_inter)
        
        for i in range(len(signals_inter)):
            access_cod = signals_inter.loc[i, 'Protein Accession code']
            access_cod = access_cod.split(';')
            organism = signals_inter.loc[i, 'Organism']
            organism = organism.split(';')
            protein_ = signals_inter.loc[i, 'Protein']
            protein_ = protein_.split(';')
            uni = signals_inter.loc[i, 'Unique Pep']
            uni = uni.split(';')
            
            idx2 = 0
            for j in range(len(access_cod)):
                prot_i = access_cod[j]
                if prot_i == protein_code:
                    idx2 = j
            
            signals_inter.loc[i, 'Protein Accession code'] = access_cod[idx2]
            signals_inter.loc[i, 'Organism'] = organism[idx2]
            signals_inter.loc[i, 'Protein'] = protein_[idx2]
            signals_inter.loc[i, 'Unique Pep'] = uni[idx2]
            
    table_comprobation(protein_code)
    table_examined = table_request(protein_code, signals_inter, choose = choose)
    try:
        sr.eliminate_table('Aquasearch_study', protein_code)
        sr.create_table_protein(db, protein_code)
        sr.insert_spectrum(db, table_examined, protein_code)
    except sqlite3.OperationalError:
        sr.create_table_protein(db, protein_code)
        sr.insert_spectrum(db, table_examined, protein_code)   
    
    
def table_request(protein_code, signals, ppm=100, db = 'Aquasearch_study', choose = 1):
    """This function completes the table belonging to a protein accession code 
       with the new signals in the found in the new sample  
       
       protein_code: string. The name of the table
       signals: Dataframe. The signals (mz and intensity) in the new sample 
       belonging to the protein of interes
       ppm: integer. Maximum error allowed for considering two signal as the same. 
       By default ppm = 100
       db: string. The name of the DB. By default: Aquasearch_study
       choose: integer.
       choose = 0: if a new signal is considered as the same as an database signal (according to ppm parameter)
                   only the signal with more relative intensity are considered 
       choose = 1: all the signals are introduced in the database              
       """
    try:
        table = sr.table_download(db, protein_code)
        table_length = len(table)
    except sqlite3.OperationalError:
        table_length = 0
        
    if table_length == 0:  # If table is empty, the signals of interest are introduced
        try:
            mz_rounded = round(signals.loc[:, 'mz'], 4)
            rel_int = round(signals.loc[:, 'relative intensity'], 2)
            table_new = pd.DataFrame({'mz': mz_rounded, 'relative intensity': rel_int, 'Unique': signals.loc[:, 'Unique Pep']})
        except KeyError:
            table_new = pd.DataFrame()
    else:
        if choose == 0:
            table = pd.DataFrame(table, columns=['mz', 'relative intensity', 'Unique'])
            new_mz = []
            new_uni = []
            new_int_r = []
        
            for mz_t, uni_t, int_re_t in zip(table.iloc[:, 0], table.iloc[:, 2], table.iloc[:, 1]):
            
                mz_d = []
                uni_d = []
                int_r_d = []
            
                for mz_s, uni_s, int_re_s in zip(signals.iloc[:, 0], signals.iloc[:, 5], signals.iloc[:, 6]):
                    
                    if mz_t > mz_s:
                        ppm_calculated = (1 - (mz_s / mz_t)) * 1000000
                    elif mz_s >= mz_t:
                        ppm_calculated = (1 - (mz_t / mz_s)) * 1000000
                    
                    if ppm_calculated <= ppm:
                        
                        if int_re_t >= int_re_s:
                            mz_d.append(mz_t)
                            uni_d.append(uni_t)
                            int_r_d.append(int_re_t)
                            
                        elif int_re_t < int_re_s:
                            mz_d.append(mz_s)
                            uni_d.append(uni_s)
                            int_r_d.append(int_re_s)
                    else:
                        mz_d.append(mz_t)
                        uni_d.append(uni_t)
                        int_r_d.append(int_re_t)
                    
                df_d = pd.DataFrame({'mz': mz_d, 'rela int': int_r_d, 'Unique': uni_d})
                df_d = df_d.drop_duplicates().reset_index().drop(['index'], axis = 1)
                
                max_int = df_d.iloc[:, 1].idxmax()
                new_mz.append(round(df_d.iloc[max_int, 0], 4))
                new_uni.append(df_d.iloc[max_int, 2])
                new_int_r.append(round(df_d.iloc[max_int, 1], 2))
             
            for mz_s, uni_s, int_re_s in zip(signals['mz'], signals['Unique Pep'], signals['relative intensity']):
                
                min_mz_ppm = mz_s-((ppm / 1000000) * mz_s)
                max_mz_ppm = (ppm / 1000000) * mz_s+mz_s

                c = False
                for mz_t in table['mz']:
                
                    if max_mz_ppm >= mz_t >= min_mz_ppm:
                        c = True
                        break

                if not c:
                    new_mz.append(round(mz_s, 4))
                    new_uni.append(uni_s) 
                    new_int_r.append(round(int_re_s, 2))
                    
            table_new = pd.DataFrame({'mz': new_mz, 'rel int': new_int_r, 'Unique': new_uni})
            table_new = table_new.drop_duplicates()
            table_new = table_new.sort_values('mz').reset_index().drop(['index'], axis = 1)
            
        elif choose == 1:
            table = pd.DataFrame(table, columns=['mz', 'relative intensity', 'Unique'])
            try:
                signals = pd.DataFrame({'mz': signals.iloc[:, 0], 'relative intensity': signals.iloc[:, 6], 'Unique': signals.iloc[:, 5]})
            except IndexError:
                pass
            table_new = table.append(signals, ignore_index= True)
            
            table_new['mz'] = round(table_new['mz'], 4)
            table_new['relative intensity'] = round(table_new['relative intensity'], 2)
            table_new = table_new.drop_duplicates(subset=['mz'])
            table_new = table_new.sort_values('mz').reset_index().drop(['index'], axis = 1)
            
            
    return table_new

def relat_intensity_calc(table_):
    """This function calculates the relative intensity of the signals belonging to a sample
       and add this information as a new column 
       
       table_: Dataframe. The result from pd_maldi_match.xml_complete.py
    """
    table_final = table_
    max_int = max(table_final.iloc[:, 1])
    list_ri = []

    for intens in table_final.iloc[:, 1]:
        ir = intens / max_int * 100
        list_ri.append(ir)
        
    table_final['relative intensity'] = list_ri
    
    return table_final

# To test the functions
if __name__ == '__main__':
    import pd_maldi_match as pdmm
    
    code = input('Uniprot code of the protein you want to search: ')
    test_pdmm = pdmm.xml_complete('test_files/mcE61_Figueres.xml',
                                  'test_files/mcE61_PD14_Figueres_Peptides.xlsx',
                                  'test_files/mcE61_PD14_Figueres_Proteins.xlsx')
    
    fill_table(code, test_pdmm, db='Aquasearch_study')
    
    

    