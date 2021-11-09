import carga_ficheros
import os
import sqlite3 as sql

def create_db(name):
    """"This fuction creates a new database in SQLite3:
        name: a string variable, it is the name we want to assign to the new DB """
    cone = sql.connect(name)
    cone.commit()
    cone.close()
    
def create_table(name, organism, protein):
    """"This fuction creates a new table in a databese:
        name: a string variable, it is the name of the DB 
        organism: a string variable, it is the name of the organism the data came from
        protein: a string vriable, it is the name of the identified protein"""
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE " + organism + "_" + protein + "(mz integer, intensity integer)")
    conn.commit()
    conn.close()
    
def df2ListOfTubles(df):
    """This function provides the transformation of the dataframe format (output of parse_zml function from 
        carga_ficheros.py) to lost of tuples for inserting the information into the table"""
    out_list = []
    for i in range(df.shape[0]):
        mz = int(df.iloc[i, 0])
        intensity = int(df.iloc[i, 1])
        tup = (mz, intensity)
        out_list.append(tup)
    return out_list

def insertSpectrum(name, organism, protein, df):
    """"This fuction complete a table in a databese:
        name: a string variable, it is the name of the DB 
        organism: a string variable, it is the name of the organism the data came from
        protein: a string vriable, it is the name of the identified protein
        df: a list of tuples, it contains the mz and intensity values of the spectrum"""
    conn = sql.connect(name)
    df_def = df2ListOfTubles(df)
    cursor = conn.cursor()
    instruction = f"INSERT INTO " + organism + "_" + protein + " VALUES (?,?)" 
    cursor.executemany(instruction,df_def)
    conn.commit()
    conn.close()


if __name__ == '__main__':
    a = input('Is this a new database? (yes/no): ')
    if a == 'yes': 
        create_db('reference_spectrums')
    
    b = input('Is this a new table? (yes/no): ')
    if b == 'yes': 
        create_table('reference_spectrums', 'pig', 'albumin')
    
        carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch/1'
        afile = 'peaklist.xml'

        path_ = os.path.join(carlos, afile)
    
        df, mz_int = carga_ficheros.parse_xml(path_)
    
        insertSpectrum('reference_spectrums', 'pig', 'albumin', df)