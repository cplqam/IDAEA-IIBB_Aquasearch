import sys

import numpy as np
import load_archives
import sqlite_requests as sr
import pandas as pd

# TODO: should be relocated on a common module and imported if needed
def dprint(*this):
    """To allow debugging with prints when using doctest"""
    that = ''.join(str(item) for item in this)
    print(that, file=sys.stderr)


def collect_peak_pairs(spec1: np.ndarray, spec2: np.ndarray, p_infor, s_infor,
                       tolerance, shift: float = 0, mz_power: float = 0.0,
                       intensity_power: float = 1.0, unique_: float = 1.5, 
                       standard: float = 1.5):
    # pylint: disable=too-many-arguments
    """Find matching pairs between two spectra.
    Args
    ----
    spec1:
        Spectrum peaks and intensities as np array.
    spec2:
        Spectrum peaks and intensities as np array.
    tolerance
        Peaks will be considered a match when <= tolerance apart.
    shift
        Shift spectra peaks by shift. The default is 0.
    mz_power:
        The power to raise mz to in the cosine function. The default is 0, in which
        case the peak intensity products will not depend on the m/z ratios.
    intensity_power:
        The power to raise intensity to in the cosine function. The default is 1.
    Returns
    -------
    matching_pairs : np array
        Array of found matching peaks.
    """

    matches = find_matches(spec1[:, 0], spec2[:, 0], p_infor, s_infor, tolerance, shift)
    try:
        dat_base = spec2[matches[:, 1].astype(int), :]
        dat_base = np.array([dat_base[:, 0], dat_base[:, 1], matches[:, 2], matches[:, 3], matches[:, 4]]).transpose()
        req = spec1[matches[:, 0].astype(int), :]
        req[:, 1] = (req[:, 1] / max(req[:, 1])) * 100
        if len(req) == 0:
            pass
    
        matching_pairs = []
        for i in enumerate(req):
            power_prod_spec1 = (req[i[0], 0] ** mz_power) * (req[i[0], 1] ** intensity_power)
            power_prod_spec2 = (dat_base[i[0], 0] ** mz_power) * (dat_base[i[0], 1] ** intensity_power)
            matching_pairs.append([req[i[0], 0], dat_base[i[0], 0], power_prod_spec1 * power_prod_spec2,
                                   dat_base[i[0], 2], dat_base[i[0], 3], dat_base[i[0], 4]])

        matching_pairs = np.array(matching_pairs)
        return matching_pairs
    except IndexError:
        pass


def find_matches(spec1_mz: np.ndarray, spec2_mz: np.ndarray, p_inf, s_inf,
                 tolerance, shift: float = 0):
    """Faster search for matching peaks.

    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from low to high m/z).

    Parameters
    ----------
    spec1_mz:
        Spectrum peak m/z values as np array. Peak mz values must be ordered.
    spec2_mz:
        Stored spectrum peak m/z values in SQLite database as np array.
        Peak mz values must be ordered.
    tolerance
        Peaks will be considered a match when <= tolerance apart.
    shift
        Shift peaks of second spectra by shift. The default is 0.
    Returns
    -------
    matches
        List containing entries of type (idx1, idx2).
    """

    p_inf = map(int, p_inf)
    p_inf = list(p_inf)

    lowest_idx = 0
    matches = np.empty((0, 5))
    for peak1_idx in range(0, spec1_mz.shape[0]):
        mz_value = spec1_mz[peak1_idx]
        low_bound = ((tolerance/1000000)-1)*-mz_value
        high_bound = -mz_value/((tolerance/1000000)-1)
        for peak2_idx in range(lowest_idx, spec2_mz.shape[0]):
            mz2 = spec2_mz[peak2_idx] + shift
            if mz2 > high_bound:
                break
            if mz2 < low_bound:
                lowest_idx = peak2_idx
            if (mz2 < high_bound) & (mz2 > low_bound):       
                if mz_value > mz2:
                    ppm = round((1-(mz2/mz_value))*1000000, 1)
                elif mz2 > mz_value:
                    ppm = round((1-(mz_value/mz2))*1000000, 1)
                else:
                    ppm = 0
                matches = np.append(matches, np.array([[int(peak1_idx), int(peak2_idx), int(p_inf[peak2_idx]),
                                                        int(s_inf[peak2_idx]), ppm]]), axis=0)
                
    matches = np.array(matches)
    # matches = matches.astype(int)
    return matches


def score_best_matches(matching_pairs: np.ndarray, mz_power: float = 0.0, intensity_power: float = 1.0, 
                       db = 'Aquasearch_study'):
    """Calculate cosine-like score by multiplying matches.

    Does require a sorted list of matching peaks (sorted by intensity product).
    """

    score = 0.0
    used_matches = int(0)
    unique_used = int(0)
    rep = np.empty((0, 7), int)
    signals_used = np.empty((0, 7), int)
    
    uni = []
    for pep in matching_pairs[:,4]:
        uni.append(sr.consulta_uniq(db, pep))
    
    ppm = matching_pairs[:, 5]
    matching_pairs = np.delete(matching_pairs, 5, axis=1)
    matching_pairs = np.insert(matching_pairs, matching_pairs.shape[1], np.asarray(uni), axis = 1)
    matching_pairs = np.insert(matching_pairs, matching_pairs.shape[1], np.asarray(ppm), axis = 1)

    for i in range(matching_pairs.shape[0]):
        mz_1 = matching_pairs[i, 0]
        mz_2 = matching_pairs[i, 1]
        sco = matching_pairs[i, 2]
          # If the peak only appears 1 time, it is added directly to the score
        if np.count_nonzero(matching_pairs[:, 0] == mz_1) == 1 and np.count_nonzero(matching_pairs[:, 1] == mz_2) == 1:
            score += matching_pairs[i, 2]
            used_matches += 1
            unique_used += matching_pairs[i, 5]
            arr = np.array([matching_pairs[i, 0], matching_pairs[i, 1], matching_pairs[i, 2], matching_pairs[i, 3],
                            matching_pairs[i, 4], matching_pairs[i, 5], matching_pairs[i, 6]]).transpose()
            signals_used = np.concatenate((signals_used, [arr]), axis=0)
          # If not, it is added to an array where the repeats will be
        else:
            rep = np.append(rep, np.array([[mz_1, mz_2, sco, matching_pairs[i, 3], matching_pairs[i, 4],
                                            matching_pairs[i, 5], matching_pairs[i, 6]]]), axis=0)
    rep_def = np.empty((0, 7), int)

    #Within the repeated
    for i in range(rep.shape[0]):
        mz_1 = rep[i, 0]
        mz_2 = rep[i, 1]
        sco = rep[i, 2]
        prot = rep[i, 3]
        pep = rep[i, 4]
        uni_ = rep[i, 5]
        ppm_ = rep[i, 6]
        if i == 0:
            # If it is the 1 element, it enters the definitive array
            rep_def = np.append(rep_def, np.array([[mz_1, mz_2, sco, prot, pep, uni_, ppm_]]), axis=0)
        else:
            # If it is not the first, it is checked that there is no equal mz in the list, of any of the 2 spectra
            # If there is not the same mz
            if np.count_nonzero(rep_def[:, 0] == mz_1) == 0 and np.count_nonzero(rep_def[:, 1] == mz_2) == 0:
                # They enter the final array
                rep_def = np.append(rep_def, np.array([[mz_1, mz_2, sco, prot, pep, uni_, ppm_]]), axis=0)
              # If there is already an equal mz in the list,
            # I compare the scores and leave the best one in that position of the array
            elif np.count_nonzero(rep_def[:, 0] == mz_1) != 0:  # Compare with first data, if the first is:
                rep1 = np.where(rep_def[:, 0] == mz_1)          # Get the index
                # If its score is better, the values are replaced to get the best score for each signal
                if sco > rep_def[rep1, 2]:
                    rep_def[rep1, 0] = mz_1
                    rep_def[rep1, 1] = mz_2
                    rep_def[rep1, 2] = sco
                    rep_def[rep1, 3] = prot
                    rep_def[rep1, 4] = pep
                    rep_def[rep1, 5] = uni_
                    rep_def[rep1, 6] = ppm_
            elif np.count_nonzero(rep_def[:, 1] == mz_2) != 0:  # The same is done for the second mz
                # If the two mz match, then it would take the best of the 3 (the two matches and the new one)
                rep2 = np.where(rep_def[:, 1] == mz_2)
                if sco > rep_def[rep2, 2]:
                    rep_def[rep2, 0] = mz_1
                    rep_def[rep2, 1] = mz_2
                    rep_def[rep2, 2] = sco
                    rep_def[rep1, 3] = prot
                    rep_def[rep1, 4] = pep
                    rep_def[rep1, 5] = uni_
                    rep_def[rep1, 6] = ppm_
                
    for i in range(rep_def.shape[0]):
        score += rep_def[i, 2]
        used_matches += 1
        unique_used += rep_def[i, 5]
        arr = np.array([rep_def[i, 0], rep_def[i, 1], rep_def[i, 2], rep_def[i, 3], 
                        rep_def[i, 4], rep_def[i, 5], rep_def[i, 6]]).transpose()
        signals_used = np.concatenate((signals_used, [arr]), axis=0)
    
    peptides = []  #Obtain the seq of the peptides used
    uni = []
    for p in signals_used[:,4]:
        peptides.append(sr.consulta_pep_seq(db, p))
        uni.append(sr.consulta_uniq(db, p))
    
    for u in range(len(uni)):
        if uni[u] == '0':
            uni[u] = 'No'
        elif uni[u] == '1':
            uni[u] = 'Yes'
        
    peptides_result = pd.DataFrame({'Peptide seq': peptides, 'unique': uni, 'ppm': signals_used[:, 6]})
    
    if score != 0:
        norm_index = signals_used[:, 0]**mz_power * signals_used[:, 2] ** intensity_power * used_matches ** (
                      (2 * unique_used)/used_matches) * ((np.sqrt(used_matches-unique_used + 1)) / 5)
        score = np.sum(norm_index) / score * 10
        
    else:
        score = 0

    return score, used_matches, unique_used, peptides_result

def request_scores(request_txt, tolerance: float, shift: float = 0,
                   mz_power: float = 0.0, intensity_power: float = 1.0, dat_b='Aquasearch_study'):
    """This function has provides a comparation and scoring of the request MALD-TOF spectrum 
       versus each protein of the database
       
       >>> r, r2 = request_scores('test_files/mcE61_Figueres.xml', tolerance = 0.05, shift =  0, mz_power = 0,intensity_power = 1.0, dat_b='Aquasearch_study')
       >>> len(r) == len(r2)
       True
       
       >>> len(r) == len(sr.get_table_names('Aquasearch_study'))
       True
       
       INPUT
       request_txt: txt file. The path of the xml file of interest
       
    """

    request = load_archives.parse_txt_request(request_txt)
    scores = {}
    # scores_unique = {}
    # peptides = {}
    
    database = load_archives.db_table_request(dat_b)
    for n in list(database.keys()):
        table_ = database[n]
        mz_intens = np.array(table_.iloc[:, [2, 3]])
        prot_inf = list(table_.iloc[:, 0])
        seq_inf = list(table_.iloc[:, 1])
        
        try:
            matches = collect_peak_pairs(request[1], mz_intens, prot_inf, seq_inf, tolerance, shift, mz_power,
                                         intensity_power)
        
            protein = sr.consulta_prot_id(dat_b, n)
            data = sr.table_request_prot_dict(dat_b, 'protein_dictionary', protein)
            
            try:
                score_r, used_m, unique_m, peptides_res = score_best_matches(matches, mz_power, intensity_power, dat_b) 
                if not data:
                    scores[protein] = [['-','-',round(score_r,2),used_m, int(unique_m)],peptides_res]
                else:
                    scores[protein] = [[data[0][0], data[0][1], round(score_r,2),used_m, 
                                        int(unique_m)],peptides_res]
            except AttributeError:
                if not data:
                    scores[protein] = [['-', '-', 0, 0, 0],[0, 0]]
                else:
                    scores[protein] = [[data[0][0], data[0][1], 0, 0, 0],[0, 0]] 
        except AttributeError:
            pass
        except ValueError:
            protein = sr.consulta_prot_id(dat_b, n)
            data = sr.table_request_prot_dict(dat_b, 'protein_dictionary', protein)
            if not data:
                scores[protein] = [['-', '-', 0, 0, 0],[0, 0]]
            else:
                scores[protein] = [[data[0][0], data[0][1], 0, 0, 0],[0, 0]]
    return scores


if __name__ == '__main__':
    
    ####Esto irá dentro de otra función
    # database = load_archives.db_table_request('Aquasearch_study')
    # table_ = database[5]
    # mz_intens = np.array(table_.iloc[:, [2, 3]])  #######Transformar
    # prot_inf = list(table_.iloc[:, 0])
    # seq_inf = list(table_.iloc[:, 1])
    
    # query1, query2 = load_archives.parse_txt_request('test_files/mcE61_Figueres.txt')
    
    # mm = find_matches(mz_intens[:, 0], mz_intens[:, 0], prot_inf, seq_inf, tolerance=100)
    # m22 = collect_peak_pairs(mz_intens, mz_intens, prot_inf, seq_inf,
    #                         tolerance = 100, shift = 0, mz_power = 0.0,
    #                         intensity_power = 1.0)
    # m32 = score_best_matches(m22, mz_power= 0.0, intensity_power = 1.0)
    scores2 = request_scores('test_program/pmf_B2.txt', tolerance= 100, shift= 0, mz_power= 0.0, intensity_power = 1.0, dat_b='Aquasearch_study')