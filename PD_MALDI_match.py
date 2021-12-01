import pandas as pd
import os
import PD_table_selection
import carga_ficheros
import numpy


def complete_table_proteins(proteins_file,n = 10):
    """This function selects and filter the most represented proteins
        
        INPUT
        proteins_file: a string vaiable, it contains the path of the file with the identified proteins
        n = an integer variable, it contains the number of proteins selected depending on their representation"""
    
    df = pd.read_excel(proteins_file)
    df = df.head(n)
    df_description = df.loc[:, 'Description'].tolist()
    name_o = []
    name_p = []
    name_p_o = []
    
    for prot in df_description:
        n_o = prot.split('=')
        n_o[1] = n_o[1].replace(' OX', '')
        n_o[0] = n_o[0].replace(' OS', '')
        n_o_join = n_o[0] + '|' + n_o[1] + ' '
        
        name_o.append(n_o[1])
        name_p.append(n_o[0])
        name_p_o.append(n_o_join)
        
    df_final = pd.DataFrame({'Accession': df.iloc[:,0], 'Organism Name': name_o, 'Protein Name':name_p})
    list_final = name_p_o
    df = df.iloc[:, 2:]
    df_final = pd.concat([df_final,df], axis=1)
    
    return df_final, list_final


def maldi_ident_join(dictionary,maldi):
    """This function joins the result of maldi and the signal identifications
    
        INPUT    
        dictionary: a dictionary variable, it contains the mz value and the protein and organism names selected
        for that mz
        maldi: an array variable, it contains the maz and intensity from MALDI spectrum"""
        
    dic_keys = list(dictionary.keys())
    mz_maldi = maldi[:,0]

    size = mz_maldi.shape[0]
    lista = ['-'] * size

    df = pd.DataFrame({'mz': maldi[:,0], 'intensity': maldi[:,1], 'Protein': lista, 'Organism': lista})

    for mz_k in dic_keys:
        p_o = dictionary[mz_k][0].split('|')
        prot = p_o[0]
        org = p_o[1]
    
        pos = numpy.where( mz_maldi == mz_k)
        df.loc[pos[0],'Protein'] = prot
        df.loc[pos[0],'Organism'] = org
    return df


def xml_complete(xml_, ident_pep, ident_prot, n_ = 10, ppm = 100, unique_ = 1):
    """This fuction assigns the protein and organism name to the MALDI spectrum signals 
    
        INPUT
        xml_: a string variable, it contains the path of the MALDI xml format file
        ident_pep: a string vaiable, it contains the path of the file with the identified peptides from PD
        ident_prot: a string vaiable, it contains the path of the file with the identified proteins from PD
        n_: an integer variable, it contains the number of proteins selected depending on their representation. 
        By default n_ = 10
        ppm: an integer variable, the maximum error allowed for considering a signal from MALDI and Orbitrap 
        as the same. By default ppm = 100
        unique_: an integer variable, if it is 1 the non-unique peptides are not included, if it is 0 the non-unique 
        peptides are not included. By default unique_ = 1"""
    
    df, mz_int = carga_ficheros.parse_xml(xml_)
    identificaciones = pd.read_excel(excel)
    identificaciones = PD_table_selection.organism_selection(identificaciones)
    df_identif, list_ident = complete_table_proteins(prot_ident,n_)   
    
    
    diccionario = {}

    for i in range(mz_int.shape[0]):
        mz_xml = mz_int[i, 0]
        lista_po = []
    
        for j in range(identificaciones.shape[0]):
            mz_ident = identificaciones.loc[j,['MH+ [Da]', 'Unique Pep']]
        
            if unique_ == 1:
                if mz_xml > mz_ident[0] and mz_ident[1] == 'Unique':
                    ppm_calculada = (1-(mz_ident[0]/mz_xml))*1000000
                elif mz_xml <= mz_ident[0] and mz_ident[1] == 'Unique':
                    ppm_calculada = (1-(mz_xml/mz_ident[0]))*1000000
    
                if ppm_calculada <= ppm:
                    app = identificaciones.loc[j,'Protein Name']+ '|' + identificaciones.loc[j,'Organism Name']
                    lista_po.append(app)
                    
            elif unique_ == 0:
                if mz_xml > mz_ident[0] :
                    ppm_calculada = (1-(mz_ident[0]/mz_xml))*1000000
                elif mz_xml <= mz_ident[0] :
                    ppm_calculada = (1-(mz_xml/mz_ident[0]))*1000000
    
                if ppm_calculada <= ppm:
                    app = identificaciones.loc[j,'Protein Name']+ '|' + identificaciones.loc[j,'Organism Name']
                    lista_po.append(app)
    
        if len(lista_po) >= 2:
            positions = []
            for num in range(len(lista_po)):
                c = []
                query = lista_po[num]
                for num2 in range(len(list_ident)):
                    answ = list_ident[num2]
        
                    if query == answ:
                        c.append(list_ident.index(query))
                    if len(c) == 0:
                        c.append(n_+1)
        
                positions.append(c[0])
        
            min_value = min(positions)    
            prot_def_inde = positions.index(min_value)
            if positions[prot_def_inde] < n_+1:
                prot_def = lista_po[prot_def_inde] 
                lista_po = prot_def
                lista_po = [lista_po]
                diccionario[mz_xml] = lista_po  #Those peptides out of the choosen top of represented proteins are deleted
        elif len(lista_po) == 1:
            diccionario[mz_xml] = lista_po
            
    result = maldi_ident_join(diccionario,mz_int)
    
    return  result
    
    

if __name__ == '__main__':
    path_ = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Datos_Ester/MALDI_Bruker_Files/MALDI_Bruker/mcE61/mcE61_Figueres/0_F3/1/1SRef/pdata/1'
    path_2 = 'C:/Users/carlos/Desktop/Carlos/Proteomics/Estudio Aquasearch/Correo de Montse/Identificaciones'
    
    afile = 'peaklist.xml'
    afile2 = 'mcE61_Figueres_01_peptides_2.xlsx'
    afile3 = 'mcE61_Figueres_01_proteins.xlsx'
    
    xml_ = os.path.join(path_,afile)
    excel = os.path.join(path_2,afile2)
    prot_ident = os.path.join(path_2,afile3)
    
    resultado = xml_complete(xml_,excel,prot_ident)
    
    
