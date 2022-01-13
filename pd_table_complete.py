import pandas as pd 
import os
import urllib.request as ur
import sqlite_requests as sr
import sqlite3


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

    group = df.loc[:, 'Protein Group Accessions'].unique()
    name = {}
    organ = {}
    exception_list = []
    
    for prot in group:
        protein = prot.split(';')
        l_n = ''
        l_o = ''
            
        for p in protein:
            data = sr.table_request_prot_dict(db, table, p)
            if len(data) == 0:
                # noinspection PyUnresolvedReferences
                try:
                    n, o = uniprot_information(p)
                    df_u = pd.DataFrame({'Accession': [p], 'Protein name': [n], 'Organism': [o]})
                    l_n = l_n + '|' + str(n)
                    l_o = l_o + '|' + str(o)
                    sr.insert_prot_code(db, table, df_u)
                except ur.HTTPError:
                    exception_list.append(p)
                    protein.remove(p)
                    if len(protein) > 1:
                        prot = ';'.join(protein)

            else:
                n = data[0][0]
                o = data[0][1]
                l_n = l_n + '|' + str(n)
                l_o = l_o + '|' + str(o)
            
        if len(prot) >= 1:            
            name[prot] = l_n[1:]
            organ[prot] = l_o[1:]

    # Remove any possible exception because the uniprot code is unavailable 
    if len(exception_list) > 0:
        exception_list = list(dict.fromkeys(exception_list))
        print('These protein accesion codes are not recognized by Uniprot and therefore, they have been eliminated from the peptide file: ' + ', '.join(exception_list))
    
        for exception_p in exception_list:
            for n in range(df.shape[0]):
                prot = df.loc[n, 'Protein Group Accessions']
                protein = prot.split(';')
            
                if exception_p in protein:
                    protein.remove(exception_p)
                    c = 0
                    if len(protein) > 1:
                        protein = ';'.join(protein)
                        df.loc[n, 'Protein Group Accessions'] = protein
                    elif len(protein) == 1:
                        df.loc[n, 'Protein Group Accessions'] = protein
                    elif len(protein) == 0:
                        df = df.drop(n-c, axis=0)
                        c += 1        
        
        name_p = table_complete(df, name)
        name_o = table_complete(df, organ)
        df['Protein Name'] = name_p
        df['Organism Name'] = name_o
        
        df.reset_index(inplace=True, drop=False) 
        df = df.drop('index', axis=1)
    return df


# To test the functions
if __name__ == '__main__':
    path_ = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Correo de Montse/Identificaciones'
    path_2 = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/IDAEA-IIBB'
    d = 'Protein_dic.xlsx'
    excel = 'mcE61_Figueres_01_peptides.xlsx'
        
    afile = os.path.join(path_, excel) 
    
    final = protein_information(afile)
