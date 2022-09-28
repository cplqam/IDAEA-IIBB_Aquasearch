###### THIS SCRIP IS A TEST WITH SAMPLES PREPARED FROM A MIX OF OTHER SAMPLES PREVIOUSLY COLLECTED
###### ALL THESE SAMPLES ARE STORED IN PROGRAM_TEST FOLDER

import score_matchms as sm

dic = {}
for f in 'B', 'D', 'F':
    for n in range(12):
        n = n+1

        file_ = 'test_program/pmf_' + f + str(n) + '.txt'

        scores, su, pep = sm.request_scores(file_, tolerance= 0.05, shift= 0, mz_power= 0.0,
                                            intensity_power = 1.0, dat_b='Aquasearch_study')

        dic['pmf_' + f + str(n)] = (scores, su, pep)
        print(file_)
