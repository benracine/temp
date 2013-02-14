from __future__ import division
import os
import shutil
import sys

import numpy as np
from numpy.random import random


def delete_old_io_files():
    try:
        os.remove('test.in')
        os.remove('test.out')
    except:
        pass


def recompile_source():
    os.system('./install.sh')


def make_input_file(year, n):
    ofile = open('test.in', 'w')
    for i in range(n):
        ofile.write('0;0;' + str(year) + '\n')
    ofile.close()


def run_shg(year):
    # cutoff_year = year + 100
    string = './lbc_smokehist_osx64.exe data/shg2p0 ' + ' ' + \
             str(int(random() * 10000)) + ' ' + \
             str(int(random() * 10000)) + ' ' + \
             str(int(random() * 10000)) + ' ' + \
             str(int(random() * 10000)) + ' ' + \
             'test.in test.out 1 0'  # -c ' + str(cutoff_year)
    print string
    os.system(string)


def test_output_file_existence():
    try:
        lines = open('test.out', 'r').readlines()
        if len(lines) > 10:
            print str(year)
            print 'nSims = ' + str(len(lines))
        else:
            print 'Run for year ' + str(year) + " didn't work"
    except:
        print 'Run for year ' + str(year) + " didn't work"


def parse_results(year):
    global sims
    people = []
    sims = [sim.split(';')[0:-1] for sim in
            open('test.out', 'r').read().split('\n')[0:-1]]
    for sim in sims:
        person = {}
        person['race'] = sim[0]
        person['sex'] = sim[1]
        person['yob'] = sim[2]
        person['init_age'] = sim[3]
        person['cess_age'] = sim[4]
        person['ocd_age'] = sim[5]
        if person['init_age'] != '-999':
            smoking_history = np.array(sim[6:])
            years = smoking_history.shape[0] / 2
            person['smoking_history'] = smoking_history.reshape(years, 2)
        people.append(person)
    return people


def comment_on_prevalence(people):
    smokers = 0
    for person in people:
        if 'smoking_history' in person.keys():
            if person['ocd_age'] > person['smoking_history'][0][0] or person['ocd_age'] == '-999':
                smokers += 1
    print 'prevalence = ' + str(smokers / len(people))


def comment_on_ages(people, key, year):
    init_ages = [int(person[key]) for person in people]
    init_ages = [d for d in init_ages if d != -999]
    local_min = min(init_ages)
    local_max = max(init_ages)
    local_median = np.median(init_ages)

    if key is 'smoking_history':
        local_min = min(init_ages)
        local_max = max(init_ages)
        local_median = np.median(init_ages)

    print key + ' = ' + str(local_min) + ', ' + str(local_max) + ', ' + str(local_median) + ' (min, max, median)'

    # Addressing minimums
    # Minimum init age is normally within some delta of 8 years old, using 6 years here
    if key is 'init_age' and abs(local_min - 8) > 7:
        print 'WARNING: weird min init age, ' + str(local_min)
    # Minimum cess age is normally within some delta of 15 years old, using 6 years here
    if key is 'cess_age' and abs(local_min - 15) > 7:
        print 'WARNING: weird min cess age, ' + str(local_min)
    # Infant mortality in the first year is normally unavoidable
    if key is 'ocd_age' and local_min is not 0:
        print 'WARNING: weird min ocd age, ' + str(local_min)

    # Addressing maximums
    # Maximum init age is normally within some delta of 8 years old, using 6 years here
    """
    if key is 'init_age' and abs(local_max - 99) > 2:
        print 'WARNING: weird max init age, ' + str(local_min)
    """
    # Maximum cess age is normally within some delta of 15 years old, using 6 years here
    if key is 'cess_age' and abs(local_max - 99) > 9 and local_max + year != 2050:
        print 'WARNING: weird max cess age, ' + str(local_max)
    # Infant mortality in the first year is normally unavoidable
    if key is 'ocd_age' and abs(local_max - 99) > 9 and local_max + year != 2050:
        print 'WARNING: weird max ocd age, ' + str(local_max)


def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]


def comment_on_smoking_histories(people, key, year):
    smoking_histories = [person['smoking_history'] for person in people if 'smoking_history' in person.keys()]
    smoking_ages = [list(np.array(smoking_history[:, 0], dtype='int')) for smoking_history in smoking_histories]
    smoking_amounts = [list(np.array(smoking_history[:, 1], dtype='float')) for smoking_history in smoking_histories]
    smoking_averages = [np.mean(np.array(smoking_amount)) for smoking_amount in smoking_amounts]
    num_switchers = len([d for d in smoking_averages if float(d) != int(d)])
    min_smoking_age = min(flatten(smoking_ages))
    max_smoking_age = max(flatten(smoking_ages))
    min_smoking_amount = min(flatten(smoking_amounts))
    average_smoking_amount = np.mean(np.array(flatten(smoking_amounts)))
    max_smoking_amount = max(flatten(smoking_amounts))
    fraction_of_switchers = num_switchers / len(smoking_histories)
    print 'Smoking ages: min = {min_smoking_age}, max = {max_smoking_age}'.format(**locals())
    print 'Smoking amounts: min = {min_smoking_amount}, max = {max_smoking_amount}, mean = {average_smoking_amount}'.format(**locals())
    print 'Fraction of switchers = {fraction_of_switchers}'.format(**locals())


if __name__ == '__main__':
    n = 20000
    #recompile_source()
    for year in range(1990, 2020, 10):
        print
        delete_old_io_files()
        make_input_file(year, n)
        run_shg(year)
        shutil.copy('test.out', 'test_' + str(year) + '.out')
        test_output_file_existence()
        people = parse_results(year)
        comment_on_prevalence(people)
        comment_on_ages(people, 'init_age', year)
        comment_on_ages(people, 'cess_age', year)
        comment_on_ages(people, 'ocd_age', year)
        comment_on_smoking_histories(people, 'ocd_age', year)
    print
    os.system("rm test* out*")
