import pandas as pd 
import os
import urllib.request as ur
import sqlite_requests as sr
import sqlite3
import time


# This function provides a information fom uniprot about a protein
def uniprot_information(accession):
    """Provides the name and organism of origin of a peptide through a Uniprot request.

    INPUT:
        accession: string. the Uniprot peptide accession code
    """

    data = ur.urlopen("https://www.uniprot.org/uniprot/%s.txt" % accession).read()
    data = str(data)
    data = data.split('\\n')
    organism = ''
    name = None
    
    for i in data:
        n = i.split('=')
        o = i.split('   ')
        
        if n[0] == 'DE   RecName: Full':
            name = n[1].strip(';')
            name = name.split('{')[0]
        elif o[0] == 'OS':
            o = str(o[1])
            organism = organism + o
            organism = organism.strip('.')
            organism = organism.split('(')[0]
        if o[0] == 'OC' and o[1].split(';')[0] == 'Eukaryote':
            organism = organism.split('(')[0]
         
    return name, organism

# To help to complete the output form Proteome Discoverer
def table_complete(df, dic_):
    """Helps to complete the output of Proteome Discoverer with the name and the organism of origin of the proteins.

    INPUT:
        df: pandas dataframe. Result of Proteome Discoverer search.
        dic_: dictionary. The protein code with information about the organism / protein name
    """
        
    group = df['Protein Group Accessions']
    new_col = []
    for i in group:
        new_col.append(dic_[i])
    return new_col

# To complete the output form Proteome Discoverer
def protein_information(path__, db='Aquasearch_study', table='protein_dictionary'):
    """Completes the table result of Proteome Discoverer search with information about proteins (protein name
    and organism) with a Uniprot request.
    
        >>> df = protein_information('test_files/mcE61_PD14_Figueres_Peptides.xlsx')
        >>> df.shape
        (1152, 18)
    
        >>> df.loc[11, 'Protein Name']
        'Alpha-amylase 1A;Alpha-amylase 2B;Pancreatic alpha-amylase'
        
        >>> df.loc[23, 'Organism Name']
        'Homo sapiens ;Homo sapiens '

    INPUT:
        path__: string. The path of the peptide output from Proteome Discoverer
        db: string.  Name of the database we are saving the datasets. Default: "Aquasearch_study"
        table: string array. Name of the table where we are saving the dictionary of eccession codes. Default: "protein_dictionary"
    """
    df = pd.read_excel(path__)

    # Try to create a database if it doesn't exist
    if not os.path.exists(db):
        sr.create_db(db)

    # Try to create a table if it doesn't exist
    try:
        sr.create_table_proteins_dic(db, table)
    except sqlite3.OperationalError:
        pass

    accessions = df.loc[:, 'Protein Group Accessions'].unique()
    protein_name = {}
    organism_name = {}
    exception_list = []

    for accession in accessions:
        acc_code = accession.split(';')
        list_prot = ''
        list_organism = ''
        long = len(acc_code)
        cont = 0
        
        while cont != long :
            acc = acc_code[cont]
            data = sr.table_request_prot_dict(db, table, acc)
            if len(data) == 0:
                # noinspection PyUnresolvedReferences
                try:
                    nam, org = uniprot_information(acc)
                    df_u = pd.DataFrame({'Accession': [acc], 'Protein name': [nam], 'Organism': [org]})
                    list_prot = list_prot + ';' + str(nam)
                    list_organism = list_organism + ';' + str(org)
                    sr.insert_prot_code(db, table, df_u)
                except ur.HTTPError:
                    exception_list.append(acc)
                    acc_code.remove(acc)
                    if len(acc_code) > 1:
                        accession = ';'.join(acc_code)
                    cont = cont - 1
                    long = long-1
                except ur.URLError: #If the host brak the conexion: sleep 10 sec and retry
                    time.sleep(10)
                    nam, org = uniprot_information(acc)
                    df_u = pd.DataFrame({'Accession': [acc], 'Protein name': [nam], 'Organism': [org]})
                    list_prot = list_prot + ';' + str(nam)
                    list_organism = list_organism + ';' + str(org)
                    sr.insert_prot_code(db, table, df_u)

            else:
                nam = data[0][0]
                org = data[0][1]
                list_prot = list_prot + ';' + str(nam)
                list_organism = list_organism + ';' + str(org) 
            cont = cont + 1 
            
        if len(accession) >= 1:            
            protein_name[accession] = list_prot[1:]
            organism_name[accession] = list_organism[1:]
            
    # Remove any possible exception because the uniprot code is unavailable 
    if len(exception_list) > 0:
        exception_list = list(dict.fromkeys(exception_list))
        print('These protein accesion codes are not recognized by Uniprot and therefore, they have been eliminated from the peptide file: ' + ', '.join(exception_list))

        for exception_p in exception_list:
            for n in range(df.shape[0]):
                acc_code = df.loc[n, 'Protein Group Accessions']
                accession = acc_code.split(';')
            
                if exception_p in accession:
                    accession.remove(exception_p)
                    c = 0
                    if len(accession) > 1:
                        accession = ';'.join(accession)
                        df.loc[n, 'Protein Group Accessions'] = accession
                    elif len(accession) == 1:
                        df.loc[n, 'Protein Group Accessions'] = accession
                    elif len(accession) == 0:
                        df = df.drop(n-c, axis=0)
                        c += 1        
        
    name_p = table_complete(df, protein_name)
    name_o = table_complete(df, organism_name)
    df['Protein Name'] = name_p
    df['Organism Name'] = name_o
        
    df.reset_index(inplace=True, drop=False) 
    df = df.drop('index', axis=1)
    return df


# To test the functions
if __name__ == '__main__':
    
    import doctest
    doctest.testmod()
    
    final = protein_information('test_files/mcE61_PD14_Figueres_Peptides.xlsx')