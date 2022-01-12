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
    

def fill_table(protein_code, maldi_complete, db='Aquasearch_study'):
    """This function creates or complete a table with the peptide signals belonging to
       the same protein in a mix sample 
       
       protein_code: string. The name of the table (protein accession code)
       maldi_complete: Dataframe. The signals from MALDI with all peptide information.
       Result from pd_maldi_match.xml_complete()
       db: string. The name of the DB. By default: Aquasearch_study
    """
    protein_codes = maldi_complete['Protein Accession code']
    protein_codes = protein_codes.tolist()

    idx = []
    for i in range(len(protein_codes)):
        if protein_codes[i] == protein_code:
            idx.append(i)
            
    signals_inter = maldi_complete.iloc[idx, :]
    signals_inter = signals_inter.reset_index().drop(['index'], axis=1)
    
    table_comprobation(protein_code)
    table_examined = table_request(protein_code, signals_inter)
    sr.eliminate_table('Aquasearch_study', protein_code)
    sr.create_table_protein(db, protein_code)
    sr.insert_spectrum(db, table_examined, protein_code)
    
    
def table_request(protein_code, signals, ppm=100, db='Aquasearch_study'):
    """This function completes the table belonging to a protein accession code 
       with the new signals in the found in the new sample  
       
       protein_code: string. The name of the table
       signals: Dataframe. The signals (mz and intensity) in the new sample 
       belonging to the protein of interes
       ppm: integer. Maximum error allowed for considering two signal as the same. 
       By default ppm = 100
       db: string. The name of the DB. By default: Aquasearch_study
       """
    try:
        table = sr.table_download(db, protein_code)
        table_length = len(table)
    except Exception:
        table_length = 0
        raise               # to catch exception name. remove when done
        
    if table_length == 0:  # If table is empty, the signals of interest are introduced
        mz_rounded = round(signals.loc[:, 'mz'], 4)
        table_new = pd.DataFrame({'mz': mz_rounded, 'intensity': signals.iloc[:, 1]})
    else:
        table = pd.DataFrame(table)
        new_mz = []
        new_int = []
        
        for i in range(len(table)):
            mz_t = table.iloc[i, 0]
            int_t = table.iloc[i, 1]
            
            mz_d = []
            int_d = []
            
            for i2 in range(len(signals)):
                mz_s = signals.iloc[i2, 0]
                int_s = signals.iloc[i2, 1]
                
                if mz_t > mz_s:
                    ppm_calculated = (1 - (mz_s / mz_t)) * 1000000
                elif mz_s >= mz_t:
                    ppm_calculated = (1 - (mz_t / mz_s)) * 1000000
                    
                if ppm_calculated <= ppm:
                    if int_t >= int_s:
                        mz_d.append(mz_t)
                        int_d.append(int_t)
                    elif int_t < int_s:
                        mz_d.append(mz_s)
                        int_d.append(int_s)
                else:
                    mz_d.append(mz_t)
                    int_d.append(int_t)
                    
            df_d = pd.DataFrame({'mz': mz_d, 'int': int_d})
            df_d = df_d.drop_duplicates().reset_index().drop(['index'], axis=1)
            max_int = df_d.iloc[:, 1].idxmax()
            new_mz.append(round(df_d.iloc[max_int, 0], 4))
            new_int.append(df_d.iloc[max_int, 1])
            
        for i in range(len(signals)):
            mz_s = signals.iloc[i, 0]
            int_s = signals.iloc[i, 1]
            
            min_mz_ppm = mz_s-((ppm/1000000)*mz_s)
            max_mz_ppm = (ppm/1000000)*mz_s+mz_s
            c = 0                      
            
            for i2 in range(len(table)):
                mz_t = table.iloc[i2, 0]
                
                if mz_t <= max_mz_ppm and mz_t >= min_mz_ppm:
                    c =+ 1
            if c == 0:
                new_mz.append(round(mz_s, 4))
                new_int.append(int_s)            
                    
        table_new = pd.DataFrame({'mz': new_mz, 'int': new_int})
        table_new = table_new.drop_duplicates()
        table_new = table_new.sort_values('mz').reset_index().drop(['index'], axis=1)
        
    return table_new

# To test the functions
if __name__ == '__main__':
    code = input('Uniprot code of the protein you want to search: ')
    fill_table(code, test_pdmm, db='Aquasearch_study')
    

    