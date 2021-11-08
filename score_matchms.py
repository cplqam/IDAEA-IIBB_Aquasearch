from typing import List 
from typing import Tuple 
import numba
import numpy as np


@numba.njit 
def collect_peak_pairs(spec1: np.ndarray, spec2: np.ndarray,
                       tolerance: float, shift: float = 0, mz_power: float = 0.0,
                       intensity_power: float = 1.0):
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
    matches = find_matches(spec1[:, 0], spec2[:, 0], tolerance, shift) 
    idx1 = [x[0] for x in matches] 
    idx2 = [x[1] for x in matches] 
    if len(idx1) == 0: 
        return None
    matching_pairs = [] 
    for i, idx in enumerate(idx1): 
        power_prod_spec1 = (spec1[idx, 0] ** mz_power) * (spec1[idx, 1] ** intensity_power) 
        power_prod_spec2 = (spec2[idx2[i], 0] ** mz_power) * (spec2[idx2[i], 1] ** intensity_power)
        matching_pairs.append([idx1[i], idx2[i], power_prod_spec1 * power_prod_spec2])
    return np.array(matching_pairs.copy())


@numba.njit
def find_matches(spec1_mz: np.ndarray, spec2_mz: np.ndarray,
                 tolerance: float, shift: float = 0) -> List[Tuple[int, int]]:
    """Faster search for matching peaks.
    Makes use of the fact that spec1 and spec2 contain ordered peak m/z (from
    low to high m/z).
    Parameters
    ----------
    spec1_mz:
        Spectrum peak m/z values as np array. Peak mz values must be ordered.
    spec2_mz:
        Spectrum peak m/z values as np array. Peak mz values must be ordered.
    tolerance
        Peaks will be considered a match when <= tolerance apart.
    shift
        Shift peaks of second spectra by shift. The default is 0.
    Returns
    -------
    matches
        List containing entries of type (idx1, idx2).
    """
    lowest_idx = 0 
    matches = []
    for peak1_idx in range(spec1_mz.shape[0]): 
        mz = spec1_mz[peak1_idx] 
        low_bound = mz - tolerance 
        high_bound = mz + tolerance
        for peak2_idx in range(lowest_idx, spec2_mz.shape[0]):  
            mz2 = spec2_mz[peak2_idx] + shift  
            if mz2 > high_bound: 
                break
            if mz2 < low_bound:  
                lowest_idx = peak2_idx
            else:  
                matches.append((peak1_idx, peak2_idx))
    return matches


@numba.njit(fastmath=True)
# this function provides some of the best scores
def score_best_matches(matching_pairs: np.ndarray, spec1: np.ndarray, 
                       spec2: np.ndarray, mz_power: float = 0.0,
                       intensity_power: float = 1.0) -> Tuple[float, int]:
    """Calculate cosine-like score by multiplying matches.
    Does require a sorted list of matching peaks (sorted by intensity product).
    """
    score = float(0.0)
    used_matches = int(0)
    used1 = set()
    used2 = set()
    rep = np.empty((0, 3), int)

    for i in range(matching_pairs.shape[0]):    
        mz_1 = matching_pairs[i, 0]
        mz_2 = matching_pairs[i, 1]
        s = matching_pairs[i, 2]
        # If the peak only appears 1 time, it is added directly to the score
        if np.count_nonzero(matching_pairs[:, 0] == mz_1) == 1 and np.count_nonzero(matching_pairs[:, 1] == mz_2) == 1:
            score += matching_pairs[i, 2] 
            used1.add(matching_pairs[i, 0])  
            used2.add(matching_pairs[i, 1])  
            used_matches += 1
        # If not, it is added to an array where the repeats will be
        else:
            rep = np.append(rep, np.array([[mz_1, mz_2, s]]), axis=0)
        
    rep_def = np.empty((0, 3), int)
    
    # Within the repeated
    for i in range(rep.shape[0]):
        mz_1 = rep[i, 0]
        mz_2 = rep[i, 1]
        s = rep[i, 2]
        if i == 0:
            # If it is the 1 element, it enters the definitive array
            rep_def = np.append(rep_def, np.array([[mz_1, mz_2, s]]), axis=0)
        else: 
            # If it is not the first, it is checked that there is no equal mz in the list, of any of the 2 spectra
            # If there is not the same mz
            if np.count_nonzero(rep_def[:, 0] == mz_1) == 0 and np.count_nonzero(rep_def[:, 1] == mz_2) == 0:
                # They enter the final array
                rep_def = np.append(rep_def, np.array([[mz_1, mz_2, s]]), axis=0)
                
            # If there is already an equal mz in the list,
            # I compare the scores and leave the best one in that position of the array
            elif np.count_nonzero(rep_def[:, 0] == mz_1) != 0:  # Compare with first data, if the first is:
                rep1 = np.where(rep_def[:, 0] == mz_1)           # Get the index
                # If its score is better, the values are replaced to get the best score for each signal
                if s > rep_def[rep1, 2]:
                    rep_def[rep1, 0] = mz_1
                    rep_def[rep1, 1] = mz_2
                    rep_def[rep1, 2] = s
            
            elif np.count_nonzero(rep_def[:, 1] == mz_2) != 0:  # The same is done for the second mz
                # If the two mz match, then it would take the best of the 3 (the two matches and the new one)
                rep2 = np.where(rep_def[:, 1] == mz_2)
                if s > rep_def[rep2, 2]:
                    rep_def[rep2, 0] = mz_1
                    rep_def[rep2, 1] = mz_2
                    rep_def[rep2, 2] = s

    for i in range(rep_def.shape[0]): 
        score += rep_def[i, 2] 
        used1.add(rep_def[i, 0])
        used2.add(rep_def[i, 1])  
        used_matches += 1         

    # Normalize score:
    spec1_power = spec1[:, 0] ** 0 * spec1[:, 1] ** 1
    spec2_power = spec2[:, 0] ** 0 * spec2[:, 1] ** 1

    score = score/(np.sum(spec1_power ** 2) ** 0.5 * np.sum(spec2_power ** 2) ** 0.5)
    return score, used_matches
