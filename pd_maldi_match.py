import pandas as pd
import pd_table_selection
import load_archives
import numpy
pd.options.mode.chained_assignment = None

def split_description(x):
    """Split uniprot protein description in individual items

    description: "Albumin OS=Gallus gallus OX=9031 GN=ALB PE=1 SV=2 - [ALBU_CHICK]"
    returns items: "Albumin", "Gallus gallus"

    """
    description = x['Description']
    items = description.split(' OX')[0]
    protein, species = items.split(" OS=")
    return protein, species


def complete_table_proteins(proteins_file, n=10):
    """Selects and filter the most represented proteins
        
        INPUT
        proteins_file: string. Path of the file with the identified proteins
        n = integer, Number of proteins selected depending on their representation
    """
    
    df_final = pd.read_excel(proteins_file)
    df_final = df_final.head(n)

    df_final[['Protein Name', 'Organism Name']] = df_final.apply(split_description, axis=1, result_type='expand')
    list_final = df_final[['Accession', 'Organism Name', 'Protein Name']]

    # TODO: Both are dataframes. Using 'list' creates confusion. Change return names in function calls
    #      for example: df_full, df_names
    return df_final, list_final


def maldi_ident_join(dictionary, maldi):
    """This function joins the result of maldi and the signal identifications
    
        INPUT    
        dictionary: a dictionary variable, it contains the mz value and the protein and organism names selected
        for that mz
        maldi: an array variable, it contains the maz and intensity from MALDI spectrum"""
        
    dic_keys = list(dictionary.keys())
    mz_maldi = maldi[:, 0]

    size = mz_maldi.shape[0]
    alist = ['-'] * size

    df = pd.DataFrame({'mz': maldi[:, 0], 'intensity': maldi[:, 1], 'Protein': alist,
                       'Organism': alist, 'Protein Accession code': alist, 'Unique Pep': alist
                       })

    for mz_k in dic_keys:
        p_o = dictionary[mz_k][0].split('|')
        prot = p_o[0]
        org = p_o[1]
        cod = p_o[2]
        pep_ = p_o[3]
    
        pos = numpy.where( mz_maldi == mz_k)
        df.loc[pos[0], 'Protein'] = prot
        df.loc[pos[0], 'Organism'] = org
        df.loc[pos[0], 'Protein Accession code'] = cod
        df.loc[pos[0], 'Unique Pep'] = pep_
    return df


def xml_complete(xml_, ident_pep, ident_prot, n_=250, ppm=100, unique_=1):
    """Assigns the protein and organism name to the MALDI spectrum signals
    
        INPUT
        xml_: string. Path of the MALDI xml format file
        ident_pep: string. Path of the file with the identified peptides from PD
        ident_prot: string. Path of the file with the identified proteins from PD
        n_: integer. Number of proteins selected depending on their representation.
            default n_ = 10
        ppm: integer. Maximum error allowed for considering a signal from MALDI and Orbitrap  as the same.
             default ppm = 100
        unique_: integer.
                 unique_ = 1, all the possibilities are considered in non-unique peptides
                 unique_ = 0, the non-unique peptides are included.
                 default unique_ = 1
    """
    
    df, mz_int = load_archives.parse_xml(xml_)
    df_ident, list_ident = complete_table_proteins(ident_prot, n_)

    dictionary = {}
    if unique_ == 1:
        identifications = pd_table_selection.organism_selection(ident_pep, sel=2)
        
        for i in range(mz_int.shape[0]):
            mz_xml = mz_int[i, 0]
            list_po = []
            for j in range(identifications.shape[0]):
                mz_ident = identifications.loc[j, 'MH+ [Da]']
                
                # Put in a list all Orbitrap signals, with their identification,
                # which are related with the MALDI signal.
                if mz_xml > mz_ident:
                    ppm_calculated = (1 - (mz_ident / mz_xml)) * 1000000
                else:
                    ppm_calculated = (1 - (mz_xml / mz_ident)) * 1000000

                if ppm_calculated <= ppm:
                    app = '|'.join([identifications.loc[j, 'Protein Name'],
                                    identifications.loc[j, 'Organism Name'],
                                    identifications.loc[j, 'Protein Group Accessions'],
                                    identifications.loc[j, 'Unique Pep']])
                    app = app.replace(';', '|').split('|')

                    options = (len(app)-1) // 3    # Last column Unique/No Unique
                    
                    if options == 1:
                        list_po.append('|'.join(app))
                    elif options >= 2:
                        for n in range(options):
                            candidate = app[n:-1:options]
                            candidate.append('No unique')
                            a = '|'.join(candidate)
                            list_po.append(a)
                            
            list_po = list(pd.unique(list_po))              
            if len(list_po) >= 2:
                p_n = []
                p_o = []
                p_c = []
                p_u = []
                
                for n in range(len(list_po)):
                    candidate = list_po[n].split('|')
                    query = candidate[2]
                    for num2 in range(len(list_ident)):
                        answ = list_ident.iloc[num2, 0]
                        if query == answ: 
                            p_n.append(candidate[0])
                            p_o.append(candidate[1])
                            p_c.append(candidate[2])
                            p_u.append('No unique')

                if len(p_n) >= 1: 
                    p_n = ';'.join(p_n)
                    p_o = ';'.join(p_o)
                    p_c = ';'.join(p_c)
                    p_u = ';'.join(p_u)
                                
                    list_po = [p_n, p_o, p_c, p_u]
                    list_po = ['|'.join(list_po)]
                    dictionary[mz_xml] = list_po  
                
            elif len(list_po) == 1:
                 
                candidate = list_po[0].split('|')
                query = candidate[2]
                for num2 in range(len(list_ident)):
                    answ = list_ident.iloc[num2, 0]
                    
                    if query == answ: 
                        list_po = [candidate[0], candidate[1], candidate[2], candidate[3]]
                        list_po = ['|'.join(list_po)]
                        dictionary[mz_xml] = list_po  
                
    elif unique_ == 0:
        identifications = pd_table_selection.organism_selection(ident_pep, sel=1)
        for i in range(mz_int.shape[0]):
            mz_xml = mz_int[i, 0]
            list_po = []
            for j in range(identifications.shape[0]):
                mz_ident = identifications.loc[j, 'MH+ [Da]']
                if mz_xml > mz_ident:
                    ppm_calculated = (1 - (mz_ident / mz_xml)) * 1000000
                else:
                    ppm_calculated = (1 - (mz_xml / mz_ident)) * 1000000

                if ppm_calculated <= ppm:
                    app = '|'.join([identifications.loc[j, 'Protein Name'],
                                    identifications.loc[j, 'Organism Name'],
                                    identifications.loc[j, 'Protein Group Accessions'],
                                    identifications.loc[j, 'Unique Pep']
                                    ])
                    list_po.append(app)

            if len(list_po) >= 2:
                positions = []
                for query in list_po:
                    c = []
                    query = query.split('|')[2]

                    for num2 in range(len(list_ident)):
                        answ = list_ident.iloc[num2, 0]
                        if query == answ:
                            c.append(list(list_ident.iloc[:, 0]).index(query))
                            
                    if len(c) == 0:
                        c.append(n_+1)
                        
                    positions.append(c[0])
                    
                idx = positions.index(min(positions))
                if positions[idx] < n_ + 1:
                    prot_def = list_po[idx]
                    list_po = prot_def
                    list_po = [list_po]
                    dictionary[mz_xml] = list_po  
                    
            elif len(list_po) == 1:
                position = []
                query = list_po[0]
                query = query.split('|')[2]

                for num2 in range(len(list_ident)):
                    answ = list_ident.iloc[num2, 0]
                    if query == answ:
                        position.append(list(list_ident.iloc[:, 0]).index(query))
                        
                if len(position) == 0:
                    position.append(n_+1)
                    
                if position[0] < n_ + 1:
                    dictionary[mz_xml] = list_po
           
    result = maldi_ident_join(dictionary, mz_int)
    
    return result
    

if __name__ == '__main__':
    
    result_1 = xml_complete('test_files/mcE61_Figueres.xml',
                            'test_files/mcE61_PD14_Figueres_Peptides.xlsx',
                            'test_files/mcE61_PD14_Figueres_Proteins.xlsx', n_=100, unique_=1)

    print(result_1)

# C:\Python38\python.exe C:/Python38/programas/aquasearch/pd_maldi_match.py
#              mz     intensity  ... Protein Accession code           Unique Pep
# 0    832.300294   9016.435667  ...                      -                    -
# 1    927.475515   5074.444521  ...                      -                    -
# 2   1071.546896  10192.549118  ...                      -                    -
# 3   1161.503095  16112.227734  ...          P06731;P80566  No unique;No unique
# 4   1181.560038  21171.262371  ...                      -                    -
# ..          ...           ...  ...                    ...                  ...
# 64  2634.326469  23260.389467  ...                      -                    -
# 65  2812.361993  10658.778184  ...                      -                    -
# 66  2990.446068  17896.442913  ...          P0DUB6;P19961  No unique;No unique
# 67  3003.349389   6511.268556  ...                      -                    -
# 68  3016.452827  30834.559000  ...                      -                    -
#
# [69 rows x 6 columns]
#
# Process finished with exit code 0
