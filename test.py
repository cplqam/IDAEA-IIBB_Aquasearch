import load_archives
import pd_table_selection 
import pd_table_complete as pdis
import pd_maldi_match as pdmm
import protein_signals as ps
import sqlite_requests as sr
import score_matchms as sm
import pandas as pd


m = input('Which module do you want to test? (1,2,3,4,t or eg): ')

# If you select the module 1 (1) 
if m == '1':
    
    test_pdis = pdis.protein_information('test_files/mcE61_PD14_Figueres_Peptides.xlsx')
    
    test_pd_sel = pd_table_selection.organism_selection('test_files/mcE61_PD14_Figueres_Peptides.xlsx', sel = 1)


# If you select the module 2 (2) 
elif m == '2':
    
    test_xml1, test_xml2 = load_archives.parse_xml('test_files/mcE61_Figueres.xml')
    
    test_pdmm = pdmm.xml_complete('test_files/mcE61_Figueres.xml',
                                  'test_files/mcE61_PD14_Figueres_Peptides.xlsx',
                                  'test_files/mcE61_PD14_Figueres_Proteins.xlsx',  
                                  n_ = 10, ppm = 100, unique_ = 1)


# If you select the module 3 (3)
elif m == '3':
    test_pdmm = pdmm.xml_complete('test_files/mcE67_Figueres.xml',
                                  'test_files/mcE67_PD14_Figueres_Peptides.xlsx',
                                  'test_files/mcE67_PD14_Figueres_Proteins.xlsx',  
                                  n_ = 10, ppm = 100, unique_ = 1)
    
    code = input('Uniprot code of the protein you want to search: ')  # Ejm: P19121 (chicken albumin)
    ps.fill_table(code, test_pdmm, db='Aquasearch_study', options = 1)
    
    table = sr.table_download('Aquasearch_study', code)
    table = pd.DataFrame(table, columns=('mz', 'intensity', 'Unique'))
    print(table)
    
    identifications = pdmm.xml_complete('test_files/mcE61_Granollers.xml',
                              'test_files/mcE61_PD14_Granollers_Peptides.xlsx',
                              'test_files/mcE61_PD14_Granollers_Proteins.xlsx', n_ = 20, unique_= 1)
    
    ps.table_union('Rodent_albumin', 'P00687','P00689', identifications)
    table = sr.table_download('Aquasearch_study', 'Rodent_albumin')
    table = pd.DataFrame(table, columns=('mz', 'intensity', 'Unique'))
    print(table)

# If you select the module 4 (4)
elif m == '4':
    res, res2 = sm.request_scores('test_files/mcE61_Mataro.xml', tolerance = 0.05, shift =  0,
                               mz_power = 0,intensity_power = 1.0, dat_b='Aquasearch_study')
# If you select all modules (t)     
elif m == 't':

    test_pdis = pdis.protein_information('test_files/mcE61_PD14_Figueres_Peptides.xlsx')
    
    test_pd_sel = pd_table_selection.organism_selection('test_files/mcE61_PD14_Figueres_Peptides.xlsx', sel = 1)
    
    test_xml1, test_xml2 = load_archives.parse_xml('test_files/mcE61_Figueres.xml')
    
    test_pdmm = pdmm.xml_complete('test_files/mcE61_Figueres.xml', 
                                  'test_files/mcE61_PD14_Figueres_Peptides.xlsx', 
                                  'test_files/mcE61_PD14_Figueres_Proteins.xlsx',  
                                  n_ = 10, ppm = 100, unique_ = 1)
    
    code = input('Uniprot code of the protein you want to search: ')
    ps.fill_table(code, test_pdmm, db='Aquasearch_study', options = 1)  # Ejm: P19121 (chicken albumin)
    
    table = sr.table_download('Aquasearch_study', code)
    table = pd.DataFrame(table, columns=('mz', 'intensity'))
    print(table)
    
    identifications = pdmm.xml_complete('test_files/mcE61_Granollers.xml',
                              'test_files/mcE61_PD14_Granollers_Peptides.xlsx',
                              'test_files/mcE61_PD14_Granollers_Proteins.xlsx', n_ = 20, unique_= 1)
    
    ps.table_union('Rodent_albumin', 'P00687','P00689', identifications)
    table = sr.table_download('Aquasearch_study', 'Rodent_albumin')
    table = pd.DataFrame(table, columns=('mz', 'intensity', 'Unique'))
    print(table)
    res, res2 = sm.request_scores('test_files/mcE61_Figueres.xml', tolerance = 0.05, shift =  0,
                               mz_power = 0,intensity_power = 1.0, dat_b='Aquasearch_study')


# If you select a complete example with the protein P19121 (chicken albumin) from 3 season in 1 WWTP (eg)
elif m == 'eg':
    maldi = ['test_files/mcE61_Figueres.xml',
             'test_files/mcE67_Figueres.xml',
             'test_files/mcE72_Figueres.xml']
    peptides = ['test_files/mcE61_PD14_Figueres_Peptides.xlsx',
                'test_files/mcE67_PD14_Figueres_Peptides.xlsx',
                'test_files/mcE72_PD14_Figueres_Peptides.xlsx']
    proteins = ['test_files/mcE61_PD14_Figueres_Proteins.xlsx',
                'test_files/mcE61_PD14_Figueres_Proteins.xlsx',
                'test_files/mcE61_PD14_Figueres_Proteins.xlsx']
    code = input('Uniprot code of the protein you want to search: ')  # Ejm: P19121 (chicken albumin)
    
    for mal, pep, pro in zip(maldi, peptides, proteins):
        
        test_pdmm = pdmm.xml_complete(mal, pep, pro, n_ = 200, ppm = 100, unique_ = 1)
        ps.fill_table(code, test_pdmm, db='Aquasearch_study', options = 1)
    
    table = sr.table_download('Aquasearch_study', code)
    table = pd.DataFrame(table, columns=('mz', 'intensity', 'relative intensity'))
    print(table)
    
    res, res2 = sm.request_scores('test_files/mcE61_Figueres.xml', tolerance = 0.05, shift =  0,
                               mz_power = 0,intensity_power = 1.0, dat_b='Aquasearch_study')
    resb, res2b = sm.request_scores('test_files/mcE67_Figueres.xml', tolerance = 0.05, shift =  0,
                               mz_power = 0,intensity_power = 1.0, dat_b='Aquasearch_study')
    resc, res2c = sm.request_scores('test_files/mcE72_Figueres.xml', tolerance = 0.05, shift =  0,
                               mz_power = 0,intensity_power = 1.0, dat_b='Aquasearch_study')