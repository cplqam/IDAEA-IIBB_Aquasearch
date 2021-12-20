import load_archives
import pd_table_selection 
import pd_table_complete as pdis
import pd_maldi_match as pdmm


m = input('Which module do you want to test? (1,2 or t): ')

# If you select the module 1 (1) 
if m == '1':
    
    test_pdis = pdis.protein_information('test_files/mcE61_Figueres_01_peptides.xlsx')
    
    test_pd_sel = pd_table_selection.organism_selection('test_files/mcE61_Figueres_01_peptides.xlsx')

# If you select the module 2 (2) 
elif m == '2':
    
    test_xml1,test_xml2 = load_archives.parse_xml('test_files/peaklist.xml')
    
    test_pdmm = pdmm.xml_complete('test_files/peaklist.xml', 'test_files/mcE61_Figueres_01_peptides_2.xlsx', 'test_files/mcE61_Figueres_01_proteins.xlsx')

# If you select all modules (t)     
elif m == 't':

    test_pdis = pdis.protein_information('test_files/mcE61_Figueres_01_peptides.xlsx')
    
    test_pd_sel = pd_table_selection.organism_selection('test_files/mcE61_Figueres_01_peptides.xlsx')
    
    test_xml1,test_xml2 = load_archives.parse_xml('test_files/peaklist.xml')
    
    test_pdmm = pdmm.xml_complete('test_files/peaklist.xml', 'test_files/mcE61_Figueres_01_peptides.xlsx', 'test_files/mcE61_Figueres_01_proteins.xlsx')
    