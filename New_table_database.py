###### THIS SCRIP IS A WORKFLOW FOR BUILDING THE REFERENCE DATABASE WITH THE SAMPLES OF THE 10
###### WASTEWATER COLLECTED AT 3 DIFFERENT TIMES AND THE STANDARS OF THE P07724, P02770, P02768,
###### P08835, P19121 AND P49065
###### ALL THESE SAMPLES ARE STORED IN test_files FOLDER


import pd_maldi_match as pdmm
import protein_signals as ps
import standard_incorporation as si


path_ = 'test_files'
path_1_stand = 'test_files/Standares/pmf_H1.txt'
path_2_stand = 'test_files/Standares/Aquasearch_Proteins_Unique.xlsx'

m = input('Do you want to enclose 2 proteines in the same table? y/n: ')

if m == 'n':

    ###### MIX   ####################################################
    codes = ['P04746','P08835','P0DUB6','P19961','P02769','P19121','P02768','P14639',
             'P35747','P07724','P00687','P02770','Q5XLE4','P00689','P49065']

    for code in codes:
        print(code)
        for folder in 'mcE61', 'mcE67', 'mcE72':

            for file in ['Banyoles', 'Besos', 'Figueres', 'Girona', 'Granollers', 'Igualada',
                         'Manresa', 'Mataro', 'Olot', 'Vic']:
                mal = path_ + '/' + folder + '_' + file + '.txt'
                pep = path_ + '/' + folder + '_PD14_' + file + '_Peptides.xlsx'
                pro = path_ + '/' + folder + '_PD14_' + file + '_Proteins.xlsx'
                code_sample = folder + '_' + file
                print(folder + '_' + file)

                test_pdmm = pdmm.txt_complete(mal, pep, pro, n_ = 200, ppm = 100, unique_ = 1)
                ps.fill_table(code, test_pdmm, code_sample, db='Aquasearch_study', options = 1)
                print(1)

    ###### STANDARDS ##########################################
    si.standard_complete('test_files/Standares/pmf_H1.txt', path_2_stand,'P07724',
                         'standard_mouse_alb_1')
    si.standard_complete('test_files/Standares/pmf_H2.txt', path_2_stand,'P07724',
                         'standard_mouse_alb_2')
    si.standard_complete('test_files/Standares/pmf_H3.txt', path_2_stand,'P07724',
                         'standard_mouse_alb_3')
    si.standard_complete('test_files/Standares/pmf_H4.txt', path_2_stand,'P02770',
                         'standard_rat_alb_1')
    si.standard_complete('test_files/Standares/pmf_H5.txt', path_2_stand,'P02770',
                         'standard_rat_alb_2')
    si.standard_complete('test_files/Standares/pmf_H6.txt', path_2_stand,'P02770',
                         'standard_rat_alb_3')
    si.standard_complete('test_files/Standares/pmf_H7.txt', path_2_stand,'P02768',
                         'standard_human_alb_1')
    si.standard_complete('test_files/Standares/pmf_H8.txt', path_2_stand,'P02768',
                         'standard_human_alb_2')
    si.standard_complete('test_files/Standares/pmf_H9.txt', path_2_stand,'P02768',
                         'standard_human_alb_3')
    si.standard_complete('test_files/Standares/pmf_H10.txt', path_2_stand,'P08835',
                         'standard_pig_alb_1')
    si.standard_complete('test_files/Standares/pmf_H11.txt', path_2_stand,'P08835',
                         'standard_pig_alb_2')
    si.standard_complete('test_files/Standares/pmf_H12.txt', path_2_stand,'P08835',
                         'standard_pig_alb_3')
    si.standard_complete('test_files/Standares/pmf_H13.txt', path_2_stand,'P19121',
                         'standard_chicken_alb_1')
    si.standard_complete('test_files/Standares/pmf_H14.txt', path_2_stand,'P19121',
                         'standard_chicken_alb_2')
    si.standard_complete('test_files/Standares/pmf_H15.txt', path_2_stand,'P19121',
                         'standard_chicken_alb_3')
    si.standard_complete('test_files/Standares/pmf_H16.txt', path_2_stand,'P02768',
                         'standard_rabbit_alb_1')
    si.standard_complete('test_files/Standares/pmf_H17.txt', path_2_stand,'P49065',
                         'standard_rabbit_alb_2')
    si.standard_complete('test_files/Standares/pmf_H18.txt', path_2_stand,'P49065',
                         'standard_rabbit_alb_3')


elif m == 'y':

####### MIX   ####################################################
    codes = list(input('Introduce the name of the 2 proteins: ').split())
    # Ejm: P07724 P02770
    new = input('Introduce a name for the table which encloses the proteins of interes: ')
    # Ejm: Murid albumin

    for folder in 'mcE61', 'mcE67', 'mcE72':
        for file in ['Banyoles', 'Besos', 'Figueres', 'Girona', 'Granollers',
                     'Igualada', 'Manresa', 'Mataro', 'Olot', 'Vic']:
            mal = path_ + '/' + folder + '_' + file + '.txt'
            pep = path_ + '/' + folder + '_PD14_' + file + '_Peptides.xlsx'
            pro = path_ + '/' + folder + '_PD14_' + file + '_Proteins.xlsx'
            code_sample = folder + '_' + file
            print(folder + '_' + file)

            test_pdmm = pdmm.txt_complete(mal, pep, pro, n_ = 200, ppm = 100, unique_ = 1)
            ps.table_union(new, codes[0],codes[1], test_pdmm, code_sample)
            print(1)

    si.standard_complete_union('test_files/Standares/pmf_H1.txt', path_2_stand, 'P07724',
                               'Murid albumin (P07724;P02770)', 'stand_murid_albu_1')
    si.standard_complete_union('test_files/Standares/pmf_H2.txt', path_2_stand, 'P07724',
                               'Murid albumin (P07724;P02770)', 'stand_murid_albu_2')
    si.standard_complete_union('test_files/Standares/pmf_H3.txt', path_2_stand, 'P07724',
                               'Murid albumin (P07724;P02770)', 'stand_murid_albu_3')
    si.standard_complete_union('test_files/Standares/pmf_H4.txt', path_2_stand, 'P02770',
                               'Murid albumin (P07724;P02770)', 'stand_murid_albu_4')
    si.standard_complete_union('test_files/Standares/pmf_H5.txt', path_2_stand, 'P02770',
                               'Murid albumin (P07724;P02770)', 'stand_murid_albu_5')
    si.standard_complete_union('test_files/Standares/pmf_H6.txt', path_2_stand, 'P02770',
                               'Murid albumin (P07724;P02770)', 'stand_murid_albu_6')
