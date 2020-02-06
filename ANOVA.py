#%%
import itertools
import functools
import ast
from datetime import datetime

import pandas as pd
from scipy import stats
import numpy as np

#%%

# function
def LIST_OF_STRING_TO_COMMA_STRING(input_string):
    return str(input_string).replace("'",'').replace('(','').replace(')','')

def RECURSIVE_COMBINE(levels, remain_factor_list, level_dict, current_conditions):
    remain_factors = list(remain_factor_list)
    conditions = []
    for c in current_conditions:
        for l in levels:
            new_condition = c + [l]
            conditions.append(new_condition)
    if len(remain_factors) > 0:
        this_factor = remain_factors.pop()
        this_levels = level_dict[this_factor]
        RECURSIVE_COMBINE(this_levels, remain_factors, level_dict, conditions)
    else:
        return conditions

def FACTORIAL_COND(factor_list, level_dict):
    factors = list(factor_list)
    factors.reverse()
    this_factor = factors.pop()
    this_levels = level_dict[this_factor]
    first_factor_conditions = [[l] for l in this_levels]
    this_factor = factors.pop()
    this_levels = level_dict[this_factor]
    all_conditions = RECURSIVE_COMBINE(this_levels, factors, level_dict, first_factor_conditions)
    return all_conditions

def LIST_TO_LIST_OF_STRINGS(input_list):
    return [str(t) for t in input_list]

def APPEND_COMBINATIONS(result_list, combination_parameters):
    factors = combination_parameters[0]
    length = combination_parameters[1]
    result_list += list(itertools.combinations(factors,length))
    return result_list

def ALL_COMBINATION(factor_list):
    combinations = functools.reduce(
        APPEND_COMBINATIONS, 
        [(factor_list,l) for l in range(2,len(factor_list)+1)], 
        list(factor_list)
        )
    return LIST_TO_LIST_OF_STRINGS(combinations)

def FILTER_BY_CRITERIA_TUPLE(df, criteria_tuple):
    factor = criteria_tuple[0]
    level = criteria_tuple[1]
    df = df[df[factor] == level]
    return df

def FILTER_DATA(df, criteria_tuples, dv):
    filtered_data = functools.reduce(
        FILTER_BY_CRITERIA_TUPLE,
        criteria_tuples,
        df.copy()
    )
    return list(filtered_data[dv])

def CALCULATE_SS(data_list, unit_n):
    cell_sum = sum(data_list)
    this_SS = sum([x**2/unit_n for x in data_list]) - cell_sum**2 / (unit_n*len(data_list))
    return this_SS, cell_sum


# import data
DATA_FILE_NAME = 'input.txt'
data = pd.read_csv(DATA_FILE_NAME, sep='\t', lineterminator='\n')


# clean up window's mess (remove "\r")
data.columns = list(data.columns[:-1]) + [data.columns[-1].rstrip('\r')]


# readout parameters
COLUMN_N = len(data.columns)
FACTOR_N = COLUMN_N - 1 # exclude the DV column
FACTOR_LIST = list(data.columns[:-1])
TOTAL_SUBJ_N = len(data)
LEVEL_DICT = {f:list(data[f].unique()) for f in FACTOR_LIST}
CONDITION_LIST = FACTORIAL_COND(FACTOR_LIST, LEVEL_DICT)
CONDITION_N = len(CONDITION_LIST)
CONDITION_SUBJ_N = TOTAL_SUBJ_N/CONDITION_N
DEPENDENT_VARIABLE = data.columns[-1]


# check format
assert FACTOR_N > 0, 'FORMAT ERROR: No IV was found.'
assert TOTAL_SUBJ_N > 0, 'FORMAT ERROR: No subject was found.'


# print parameters
parameters_string = 'Between-group ANOVA Factors:\n'
for f in FACTOR_LIST:
    parameters_string += f + ': ' + LIST_OF_STRING_TO_COMMA_STRING(LEVEL_DICT[f]) + '\n'
parameters_string += 'DV: ' + DEPENDENT_VARIABLE
print(parameters_string)


# construct ANOVA table
results = pd.DataFrame(columns=['Source','SS','df','MS','F','p','etap^2'])

effect_list = ALL_COMBINATION(FACTOR_LIST)

source_list = ['Between']
source_list += effect_list
source_list += ['Error', 'Total']

results['Source'] = source_list
results = results.set_index('Source')


# calculate Error SS, df, and MS
sum_of_each_cell = [] # for between SS
error_SS = 0
error_df = 0
for c in CONDITION_LIST:
    assert len(FACTOR_LIST) == len(c), (
        'PROGRAM ERROR: '+
        'condition name is incorrect '+
        '(inconsistent with the number of factors specified)')
    this_data = FILTER_DATA(data, zip(FACTOR_LIST, c), DEPENDENT_VARIABLE)
    cell_SS, sum_of_this_cell = CALCULATE_SS(this_data, 1)
    error_SS += cell_SS
    error_df += len(this_data) - 1
    sum_of_each_cell.append(sum_of_this_cell) # for between SS

results.loc['Error','SS'] = error_SS
results.loc['Error','df'] = error_df
results.loc['Error','MS'] = error_SS/error_df


# calculate between SS and df
between_SS, sum_of_all_data = CALCULATE_SS(sum_of_each_cell, CONDITION_SUBJ_N)
between_df = CONDITION_N - 1
results.loc['Between','SS'] = between_SS
results.loc['Between','df'] = between_df


# calculate SS, df, and MS for main effects
subtraction_term = sum_of_all_data**2 / TOTAL_SUBJ_N
for f in FACTOR_LIST:
    this_levels = LEVEL_DICT[f]
    cell_size = TOTAL_SUBJ_N / len(this_levels)
    sum_of_cell_list = [sum(FILTER_DATA(data, [(f,l)], DEPENDENT_VARIABLE)) for l in this_levels]
    this_SS, _ = CALCULATE_SS(sum_of_cell_list, cell_size)
    this_SS -= subtraction_term
    this_df = len(this_levels) - 1
    results.loc[f,'SS'] = this_SS
    results.loc[f,'df'] = this_df
    results.loc[f,'MS'] = this_SS / this_df


# calculate SS, df, and MS for interaction effects
for s in effect_list[len(FACTOR_LIST):]:
    source_tuple = ast.literal_eval(s)
    all_dfs = [len(LEVEL_DICT[f])-1 for f in source_tuple]
    this_df = np.prod(all_dfs)

    this_SS = 0
    partition_list = ALL_COMBINATION(ast.literal_eval(s))[:-1]
    partition_SS_list = [results.loc[i,'SS'] for i in partition_list]
    partition_df_list = [results.loc[i,'df'] for i in partition_list]
    total_SS_subtraction_term = sum(partition_SS_list)

    cell_list = FACTORIAL_COND(source_tuple, LEVEL_DICT)
    for c in cell_list:
        this_data = FILTER_DATA(data, list(zip(source_tuple,c)), DEPENDENT_VARIABLE)
        cell_sum = sum(this_data)
        cell_count = len(this_data)
        this_SS += cell_sum**2/cell_count
    this_SS -= subtraction_term
    this_SS -= total_SS_subtraction_term
    results.loc[s,'SS'] = this_SS
    results.loc[s,'df'] = this_df
    results.loc[s,'MS'] = this_SS / this_df
    if s == str(tuple(FACTOR_LIST)):
        assert (
            this_SS == between_SS - sum(list(results['SS'].dropna())[1:-1]), 
            'PROGRAM ERROR: SS do not sum to total SS'
            )
        assert (
            this_df == between_df - sum(list(results['df'].dropna())[1:-1]),
            'PROGRAM ERROR: df do not sum to total df'
            )



# calculate total SS and df
results.loc['Total','SS'] = results.loc['Between','SS'] + results.loc['Error','SS']
results.loc['Total','df'] = results.loc['Between','df'] + results.loc['Error','df']


# calculate Fs, ps, and partial eta squared
for s in effect_list:
    results.loc[s,'F'] = results.loc[s,'MS'] / results.loc['Error','MS']
    results.loc[s,'p'] = 1 - stats.f.cdf(results.loc[s,'F'], results.loc[s,'df'], results.loc['Error','df'])
    results.loc[s,'etap^2'] = results.loc[s,'SS'] / (results.loc[s,'SS']+results.loc['Error','SS'])


# clean up ANOVA table
results.drop(['Between'], axis=0, inplace=True)

index_list = list(results.index)
results.index = pd.Index(
    index_list[:len(FACTOR_LIST)] + 
    [LIST_OF_STRING_TO_COMMA_STRING(s).replace(', ',' x ') for s in index_list[len(FACTOR_LIST):-2]] + 
    index_list[-2:]
    )


# print results
index_list = list(results.index)
print('\n')
for e in index_list[:-2]:
    print (e+': p = '+format(results.loc[e,'p'],'0.3f')[1:])


# save results
results.to_csv('output.txt', header=True, index=True, sep='\t')


# record parameters
with open('output.txt','a') as results_file:
    results_file.write('\n'+parameters_string)


# %%
