import carga_ficheros
import os
import sqlite3 as sql

def create_db(name):
    """"This fuction creates a new database in SQLite3:
        name: a string variable, it is the name we want to assign to the new DB """
    cone = sql.connect(name)
    cone.commit()
    cone.close()
    
#To create a protein diccionary
def create_table_proteins_dic(name, table_n):
    """"This fuction creates a new table in a databese:
        name: a string variable, it is the name of the DB 
        table_n, a string variable, it is the name of the diccionary"""
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("CREATE TABLE " + table_n + "(Accession text, Protein name text, Organism text)")
    conn.commit()
    conn.close()

def df2ListOfTubles_proteins_dic(df):
    """This function provides the transformation of the dataframe format (output of parse_zml function from 
        carga_ficheros.py) to lost of tuples for inserting the information into the table"""
    out_list = []
    for i in range(df.shape[0]):
        acces = df.iloc[i, 0]
        prot_n= df.iloc[i, 1]
        org = df.iloc[i, 2]
        tup = (acces, prot_n, org)
        out_list.append(tup)
    return out_list

def insertProtCode(name, table_n, df):
    """"This fuction complete a table in a databese:
        name: a string variable, it is the name of the DB 
        table_n, a string variable, it is the name of the diccionary
        df: a list of tuples, it contains the accession code, protein name and organism"""
    conn = sql.connect(name)
    df_def = df2ListOfTubles_proteins_dic(df)
    cursor = conn.cursor()
    instruction = f"INSERT INTO " + table_n + " VALUES (?,?,?)" 
    cursor.executemany(instruction,df_def)
    conn.commit()
    conn.close()
    
def table_request_protDic(name, table_n, code):
    """"This fuction provides to do a request to get information from the protein diccionary:
        name: a string variable, it is the name of the DB 
        table_n, a string variable, it is the name of the diccionary
        code: a string variable, it is the uniprot code of the protein"""
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("SELECT Protein name, Organism FROM " + table_n + "  WHERE Accession = '" + code + "'")
    datos = cursor.fetchall()
    conn.close()
    return datos


#To introduce a spectrum of reference
def create_table_organism(name, organism, protein):
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
    
    
def table_request_prueba(name, table_n, code):
    """"This fuction provides to do a request to get information from the protein diccionary:
        name: a string variable, it is the name of the DB 
        table_n, a string variable, it is the name of the diccionary
        code: a string variable, it is the uniprot code of the protein"""
    conn = sql.connect(name)
    cursor = conn.cursor()
    cursor.execute("SELECT intensity FROM " + table_n + "  WHERE mz = " + code + "")
    datos = cursor.fetchall()
    conn.close()
    return datos

#To test the functions
# if __name__ == '__main__':
#     a = input('Is this a new database? (yes/no): ')
#     if a == 'yes': 
#         create_db('reference_spectrums')
        
#     b = input('Is this a new table? (yes/no): ')
#     if b == 'yes': 
#         create_table('reference_spectrums', 'pig', 'albumin')
    
#         carlos = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Prueba/AquaSearch/1'
#         afile = 'peaklist.xml'

#         path_ = os.path.join(carlos, afile)
    
#         df, mz_int = carga_ficheros.parse_xml(path_)
    
#     insertSpectrum('reference_spectrums', 'pig', 'albumin', df)

    