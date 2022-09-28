import os
import pandas as pd
import sqlite_requests as sr
import numpy as np
# from statistics import mode
from collections import Counter


def table_proteins_check(protein_code, db='Aquasearch_study'):
    """Checks if the database and a table with the protein accession code name exist.
       If any of them do not exist, it is created

       protein_code: string. The name of the table (protein accession code)
       db: string. The name of the DB. By default: Aquasearch_study
       """

    # Try to create a database if it doesn't exist
    if not os.path.exists(db):
        sr.create_db(db)

    # Try to create a table if it doesn't exist or complete the table if it exists and the protein is new
    sr.protein_codes_table(db)
    table = sr.table_download(db, 'Protein_codes')
    table = pd.DataFrame(table)
    
    try:
        codes = list(table.iloc[:, 1])
        if protein_code not in codes:
            id_prot = len(codes)+1
            sr.insert_protein_code(db, protein_code, id_prot)
        else:
            id_prot = codes.index(protein_code) + 1
    except IndexError:
        id_prot = 1
        sr.insert_protein_code(db, protein_code, id_prot)
        
    return id_prot
    
def table_sequence_check(sequence_prot, uniq, m_mass, db='Aquasearch_study'):
    
    sr.table_sequences(db)
    table = sr.table_download(db, 'Peptide_sequences')
    table = pd.DataFrame(table)
    
    try:
        sequences = list(table.iloc[:, 1])
        
        if sequence_prot not in sequences:
            id_sequence = len(sequences)+1
            sr.new_sequence(db, sequence_prot, id_sequence, uniq, m_mass)
        else:
            id_sequence = sequences.index(sequence_prot) + 1
    except IndexError:
        id_sequence = 1
        sr.new_sequence(db, sequence_prot, id_sequence, uniq, m_mass)
        
    return id_sequence

def table_spectrum_check(signals, id_, db = 'Aquasearch_study'):

    sr.spectrums_table(db)
    sr.table_sequences(db)
    table_spect = sr.table_download(db, 'Spectrums_table')  
    
    table_spect = pd.DataFrame(table_spect, columns= ('mz', 'relative_intensity', 
                                                      'standard_signal', 'protein', 
                                                      'sample', 'sequence'))
    try:
        mz_rounded = round(signals.loc[:, 'mz'], 4)
        rel_int = round(signals.loc[:, 'relative_intensity'], 2)
        id_protein = [str(id_)]*signals.shape[0]
    
        table_request= pd.DataFrame({'mz': mz_rounded,
                                     'relative_intensity': rel_int,
                                     'standard_signal': signals.loc[:, 'Standard signal'],
                                     'protein': id_protein,
                                     'sample': signals.loc[:, 'Sample'],
                                     'sequence': signals.loc[:, 'Peptide seq']})
    except KeyError:
        table_request = pd.DataFrame()
    
    table_new = pd.concat([table_spect, table_request], ignore_index=True)
    table_new['mz'] = round(table_new['mz'], 4)
    table_new['relative_intensity'] = round(table_new['relative_intensity'], 2)
    table_new = table_new.drop_duplicates(subset=['mz', 'protein', 'sample', 'sequence'])
    table_new = table_new.sort_values('mz').reset_index().drop(['index'], axis=1)
    
    return table_new


def peptide_reevaluation(db='Aquasearch_study'):
    
    table = sr.table_download(db, 'Spectrums_table')
    table = pd.DataFrame(table, columns= ('mz', 'relative_intensity', 
                                          'standard_signal', 'protein', 
                                          'sample', 'sequence'))
    no_unique = []
    n_peptides = np.unique(table['sequence'])
    for pep in n_peptides:
      peptides = table[table['sequence'] == pep]
      proteins = list(np.unique(peptides['protein']))
      
      if len(proteins) > 1:
          no_unique.append(pep-1)
      
    table_pep = sr.table_download(db, 'Peptide_sequences')
    table_pep = pd.DataFrame(table_pep, columns= ('id', 'peptide_sequence', 
                                              'uniq', 'molecular_mass'))
    table_pep.loc[no_unique, 'uniq'] = '0'
    
    sr.eliminate_table(db, 'Peptide_sequences')
    sr.table_sequences(db)
    for id_sequence, sequence_prot, uniq, m_mass in zip(table_pep['id'], table_pep['peptide_sequence'], table_pep['uniq'], table_pep['molecular_mass']):
        sr.new_sequence(db, sequence_prot, id_sequence, uniq, m_mass)
      


def quality_filter(db='Aquasearch_study'):
    sr.spectrums_table_filtered(db)
    table = sr.table_download(db, 'Spectrums_table')
    table = pd.DataFrame(table, columns= ('mz', 'relative_intensity', 
                                          'standard_signal', 'protein', 
                                          'sample', 'sequence'))
    table_mix = table[table['standard_signal'] == '0']
    table_stand = table[table['standard_signal'] == '1']

    table_mix['protein'] = list(map(int,(table_mix['protein'])))

    spectrums_filtered = pd.DataFrame(columns= ('mz', 'relative_intensity', 
                                                'standard_signal', 'protein', 
                                                'sample', 'sequence'))

    n_proteins = max(np.unique(table_mix['protein']))

    for protein in range(n_proteins):
        protein = protein+1
        
        signals_protein = table_mix[table_mix['protein'] == protein]
        if signals_protein.shape[0] > 0:
            most_signals = Counter(signals_protein['sample'])
            most_signals = dict(most_signals)
        
            signals_mean = np.mean(list(most_signals.values()))
            signals_std = np.std(list(most_signals.values()))
            threshold = round(signals_mean - signals_std*0.75, 0)
        
            list_samples = np.unique(list(signals_protein['sample']))
        
            for samp in list_samples:
                each_sample = signals_protein[signals_protein['sample'] == samp]
                
                if threshold != 1:
                    if (len(each_sample) >= threshold) & (len(each_sample) > 1): #or (each_sample.iloc[0,2] == '1'):
                        spectrums_filtered = pd.concat([spectrums_filtered, each_sample], ignore_index=True)
                    else:
                        pass
                else:
                    if (len(each_sample) >= threshold):
                        spectrums_filtered = pd.concat([spectrums_filtered, each_sample], ignore_index=True)
                    else:
                        pass


    spectrums_filtered = pd.concat([spectrums_filtered, table_stand], ignore_index=True)
    spectrums_filtered = spectrums_filtered.sort_values(by=['mz'], ascending=True)


    sr.eliminate_table(db, 'Spectrums_table_filtered')
    sr.spectrums_table_filtered(db)
    sr.insert_new_spectrum_filtered(db, spectrums_filtered)
    
def reference_spectrum(db='Aquasearch_study'):
    
    sr.table_quant_inf(db)
    table = sr.table_download(db, 'Spectrums_table_filtered')
    table = pd.DataFrame(table, columns= ('mz', 'relative_intensity', 
                                          'standard_signal', 'protein', 
                                          'sample', 'sequence'))

    table['protein'] = list(map(int,(table['protein'])))

    quantitative_inf = pd.DataFrame(columns= ('protein', 'sequence', 'mz_aver', 'ir_aver'))

    protein_peptide = table[['protein', 'sequence']]
    protein_peptide = protein_peptide.drop_duplicates()

    for prot, pep in zip(protein_peptide['protein'], protein_peptide['sequence']):
        
        
        list_signals = table[table['protein'] == prot] 
        list_signals = list_signals[list_signals['sequence'] == pep]
        
        mean_mz = np.average(list_signals['mz'])
        mean_ri = np.average(list_signals['relative_intensity'])
        std_ri = np.std(list_signals['relative_intensity'])

        error_ri = []
        
        for mz, ri in zip(list_signals['mz'], list_signals['relative_intensity']):
            error_ri.append(np.sqrt((ri-mean_ri)**2))
            
        list_signals['error ri'] = error_ri
        error_ri = list(map(int,error_ri < std_ri))
        condition = 1
        
        while condition == 1:
            if list_signals.shape[0] >= 3 & sum(error_ri) >= 2:
                errores_signals = list(list_signals['error ri'])
                
                idx_ = list(map(int, errores_signals >= std_ri))
                idx_2 = idx_.index(1)
                
                if sum(idx_) >= 1:
                
                    list_signals = list_signals.drop(list_signals.index[[idx_2]])
                
                    mean_mz = np.average(list_signals['mz'])
                    mean_ri = np.average(list_signals['relative_intensity'])
                    std_ri = np.std(list_signals['relative_intensity'])
                
                    error_ri = []
                
                    for mz, ri in zip(list_signals['mz'], list_signals['relative_intensity']):
                        error_ri.append(np.sqrt((ri-mean_ri)**2))
                    
                    list_signals['error ri'] = error_ri
                    error_ri = list(map(int,error_ri < std_ri))
                    
                else:
                    condition = 0
                    
            else:
                condition = 0
        
        q_inf = pd.DataFrame({'protein': str(prot), 'sequence': str(pep), 'mz_aver': [mean_mz], 'ir_aver': [mean_ri]})
        
        quantitative_inf = pd.concat([quantitative_inf, q_inf])
        
    sr.eliminate_table(db, 'Quantitative_information') 
    sr.table_quant_inf(db)   
    sr.new_quiant_inf(db, quantitative_inf) 

def relat_intensity_calc(table_):
    """Calculates the relative intensity of the signals belonging to a sample
       and adds this information as a new column

       table_: Dataframe. The result from pd_maldi_match.txt_complete.py
    """
    table_final = table_
    max_int = max(table_final.iloc[:, 1])

    list_ri = [intens / max_int * 100 for intens in table_final.iloc[:, 1]]

    table_final['relative_intensity'] = list_ri

    return table_final
    

def fill_table(protein_code, maldi_complete, sample_name, db='Aquasearch_study', options=1):
    """Creates or completes a table with the peptide signals belonging to protein <protein_code> in a protein mix.

       protein_code: string. The name of the table (protein accession code)
       maldi_complete: Dataframe. The signals from MALDI with all peptide information.
       sample_name: string. The name of the sample from which the signals were obtained
       Result from pd_maldi_match.txt_complete()
       db: string. The name of the DB. By default: Aquasearch_study

       options: integer
       options = 0: there is only 1 peptide option for each signal
                    (pd_maldi_match.txt_complete(unique_ = 0))
       options = 1: each signal can have 1 or more peptide signals associated
                    (pd_maldi_match.txt_complete(unique_ = 1))
    """
    signals_inter = None      # signals of interest ?
    accessions = maldi_complete['Protein Accession code']
    accessions = accessions.tolist()

    if options == 0:
        idx = []
        for i, accession in enumerate(accessions):
            if accession == protein_code:
                idx.append(i)

        signals_inter = maldi_complete.iloc[idx, :]
        signals_inter = signals_inter.reset_index().drop(['index'], axis=1)
        if signals_inter.shape[0] > 0:
            signals_inter = relat_intensity_calc(signals_inter)

    elif options == 1:
        idx = []
        for i, accession in enumerate(accessions):
            acc_code = accession.split(';')

            for acc in acc_code:
                if acc == protein_code:
                    idx.append(i)

        signals_inter = maldi_complete.iloc[idx, :]
        signals_inter = signals_inter.reset_index().drop(['index'], axis=1)
        if signals_inter.shape[0] > 0:
            signals_inter = relat_intensity_calc(signals_inter)

            for i in range(signals_inter.shape[0]):
                access_cod = signals_inter.loc[i, 'Protein Accession code']
                access_cod = access_cod.split(';')
                organism = signals_inter.loc[i, 'Organism Name']
                organism = organism.split(';')
                protein_ = signals_inter.loc[i, 'Protein']
                protein_ = protein_.split(';')
                uni = signals_inter.loc[i, 'Unique Pep']
                uni = uni.split(';')
                stand_sig = signals_inter.loc[i, 'Standard signal']
                stand_sig = stand_sig.split(';')
                sequence = signals_inter.loc[i, 'Peptide seq']
                sequence = sequence.split(';')
                m_mass = signals_inter.loc[i, 'Peptide mass']
                m_mass = m_mass.split(';')

                idx2 = 0
                for j, access_j in enumerate(access_cod):
                    if access_j == protein_code:
                        idx2 = j

                signals_inter.loc[i, 'Protein Accession code'] = access_cod[idx2]
                signals_inter.loc[i, 'Organism Name'] = organism[idx2]
                signals_inter.loc[i, 'Protein'] = protein_[idx2]
                signals_inter.loc[i, 'Unique Pep'] = uni[idx2]
                signals_inter.loc[i, 'Standard signal'] = stand_sig[idx2]
                signals_inter.loc[i, 'Peptide seq'] = sequence[idx2]
                signals_inter.loc[i, 'Peptide mass'] = m_mass[idx2]
    
    
    signals_inter.loc[:, 'Sample'] = [sample_name] * signals_inter.shape[0]           
    id_p = table_proteins_check(protein_code)
    
    code_seq = []
    for seq, uniq, mass in zip(signals_inter.loc[:, 'Peptide seq'], signals_inter.loc[:, 'Unique Pep'], signals_inter.loc[:, 'Peptide mass']):
        id_s = table_sequence_check(seq, uniq, mass)
        code_seq.append(id_s)
    
    signals_inter.loc[:, 'Peptide seq'] = code_seq

    table_examined = table_spectrum_check(signals_inter, id_p)
    sr.eliminate_table(db, 'Spectrums_table')
    sr.spectrums_table(db)
    sr.insert_new_spectrum(db, table_examined)
    
    quality_filter(db)
    reference_spectrum(db)
    peptide_reevaluation(db)  


def table_union(new, old1, old2, signals, sample_name, db='Aquasearch_study'):
    """Selects the peptide signals of 2 different proteins in the same group

       old1: string. The name of one protein you want to combine
       old2: string. The name of the other protein you want to combine
       signals: Dataframe. The result from pd_maldi_match.txt_complete.py
       db: string. The name of the database
    """
    new = new + ' (' + old1 + ';' + old2 + ')'
    protein_codes = signals['Protein Accession code']
    protein_codes = protein_codes.tolist()
    mz = []
    intens = []
    uni = []
    stand = []
    seq = []
    m_mass = []
    for i, codes in enumerate(protein_codes):
        signal_code = codes.split(';')
        num_pep = len(signal_code)

        # possible counter values: [0, 1, 2]
        counter = sum([old1 in signal_code, old2 in signal_code])
        if counter == 0:
            continue
        
        elif counter == 1:
            target_signal = signals.iloc[i, :]
            codes_target = target_signal['Protein Accession code'].split(';')
            try:
                t1 = codes_target.index(old1)
                c = 1
            except ValueError:
                pass
            try:
                t2 = codes_target.index(old2)
                c = 2
            except ValueError:
                pass
            
            if c == 1:
                mz.append(signals.iloc[i,0])
                intens.append(signals.iloc[i,0])
                stand.append('0')
                seq.append(target_signal[7].split(';')[t1])
                m_mass.append(target_signal[8].split(';')[t1])
                if counter == num_pep:
                    uni.append('1')
                else:
                    uni.append('0')
            elif c == 2:
                mz.append(signals.iloc[i,0])
                intens.append(signals.iloc[i,0])
                stand.append('0')
                seq.append(target_signal[7].split(';')[t2])
                m_mass.append(target_signal[8].split(';')[t2])
                if counter == num_pep:
                    uni.append('1')
                else:
                    uni.append('0')       
            
        else:
            
            target_signal = signals.iloc[i, :]
            codes_target = target_signal['Protein Accession code'].split(';')
            t1 = codes_target.index(old1)
            t2 = codes_target.index(old2)
            
            if target_signal[7].split(';')[t1] == target_signal[7].split(';')[t2]:
                mz.append(signals.iloc[i,0])
                intens.append(signals.iloc[i,0])
                stand.append('0')
                seq.append(target_signal[7].split(';')[t1])
                m_mass.append(target_signal[8].split(';')[t1])
                if counter == num_pep:
                    uni.append('1')
                else:
                    uni.append('0')
                
            else:
                mz.append(signals.iloc[i,0])
                intens.append(signals.iloc[i,0])
                stand.append('0')
                seq.append(target_signal[7].split(';')[t1])
                m_mass.append(target_signal[8].split(';')[t1])
                if counter == num_pep:
                    uni.append('1')
                else:
                    uni.append('0')
                
                mz.append(signals.iloc[i,0])
                intens.append(signals.iloc[i,0])
                stand.append('0')
                seq.append(target_signal[7].split(';')[t2])
                m_mass.append(target_signal[8].split(';')[t2])
                if counter == num_pep:
                    uni.append('1')
                else:
                    uni.append('0')

    if len(mz) > 0:
        df_intens = pd.DataFrame({'mz': mz, 'int': intens})
        rel_intens = relat_intensity_calc(df_intens)

        rel_intens_rounded = round(rel_intens.iloc[:, 2], 2)
        mz_rounded = round(df_intens.iloc[:, 0], 4)
        signals_interest = pd.DataFrame({'mz': mz_rounded, 'relative_intensity': rel_intens_rounded,
                                          'Unique Pep': uni,'Standard signal': stand, 'Peptide seq': seq,
                                          'Peptide mass': m_mass})

        id_p = table_proteins_check(new)    
        code_seq = []
        for seq, uniq, mass in zip(signals_interest.loc[:, 'Peptide seq'], signals_interest.loc[:, 'Unique Pep'], signals_interest.loc[:, 'Peptide mass']):
            id_s = table_sequence_check(seq, uniq, mass)
            code_seq.append(id_s)

        signals_interest.loc[:, 'Peptide seq'] = code_seq
        signals_interest['Sample'] = [sample_name] * signals_interest.shape[0]


        table_examined = table_spectrum_check(signals_interest, id_p)
        sr.eliminate_table(db, 'Spectrums_table')
        sr.spectrums_table(db)
        sr.insert_new_spectrum(db, table_examined)
        
        quality_filter(db)
        reference_spectrum(db)
        peptide_reevaluation(db)


# To test the functions
if __name__ == '__main__':
    import pd_maldi_match as pdmm

    code = 'P02770'
    test_pdmm = pdmm.txt_complete('test_files/mcE72_Granollers.txt',
                                  'test_files/mcE72_PD14_Granollers_Peptides.xlsx',
                                  'test_files/mcE72_PD14_Granollers_Proteins.xlsx')

    fill_table(code, test_pdmm, sample_name='mcE61_Banyoles', db='Aquasearch_study')
    #table_union('Murid albumin','P02770', 'P07724',  test_pdmm, sample_name='mcE61_Figueres')