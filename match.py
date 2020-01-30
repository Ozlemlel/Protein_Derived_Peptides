#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 22:05:05 2019

@author: Steve
"""

import pandas as pd

from Bio import pairwise2

import matplotlib.pyplot as plt

# This interface is used for the filter final percentage function when passing in comparison operators
# Changes(Jan 28 2020): Remove filter functions for TSS >= 80 and <= 0 -> generate 2 graphs instead of four.
# Change (Jan 29 2020): 
# 1. Filter the human dataset to TSS above 80% and below 20% DONE
# 2. Filter the target dataset to TSS above 80% and below 20% DONE
# 3. Match the datasets accordingly DONE
import operator

import os

global final_percentage_score_eighty_above_ninty
global final_percentage_score_zero_above_eighty
global final_percentage_score_eighty_below_twenty
global final_percentage_score_zero_below_twenty
    
final_percentage_score_eighty_above_ninty = pd.DataFrame()
final_percentage_score_zero_above_eighty = pd.DataFrame()
final_percentage_score_eighty_below_twenty = pd.DataFrame()
final_percentage_score_zero_below_twenty = pd.DataFrame()

def get_truth(inp, relate, cut):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '=': operator.eq}
    return ops[relate](inp, cut)

def compare_data_cross_species(ref_file_name, target_file_name, matrix):    
    for filename in os.listdir(target_file_name):
        print('processing: ' + filename)
        compare_data(ref_file_name, './' + target_file_name + '/' + filename, matrix)
    # Change the second and third argument for size 
    plot_chart(final_percentage_score_eighty_above_ninty, 50, 10, 'total_above_ninty.png')
    # plot_chart(final_percentage_score_zero_above_eighty, 50, 10, 'below_zero_above_eighty_precent.png')
    plot_chart(final_percentage_score_eighty_below_twenty, 50, 10, 'total_below_twenty.png')
    # plot_chart(final_percentage_score_zero_below_twenty, 50, 10, 'below_zero_below_twenty_precent.png')
    print("DONE")
        

# main function of the program. Takes 2 datasets as input and another similarity matrix as string input.
# Sequences need to be the same, otherwise program fails
# Similarity matrix is prefered to be symmetric
# This is the only function necessary to call to get the output
# The images will be saved in the current directory where the program is located as png
# cheng the second and third argument of plot_chart to change the size of the output image
def compare_data(reference, target, matrix_input):
    global ref_eighty
    global ref_zero
    global tar_eighty
    global tar_zero
    global match_table_eighty
    global match_table_zero
    global matrix
    global min_max_eighty
    global min_max_zero
    global final_percentage_score_eighty_above_ninty
    global final_percentage_score_zero_above_eighty
    global final_percentage_score_eighty_below_twenty
    global final_percentage_score_zero_below_twenty
    
    
    ref = pd.read_csv(reference)
    tar = pd.read_csv(target)
    
    matrix_unready = pd.read_csv(matrix_input)
    matrix = prep_matrix(matrix_unready)
    
    # Finds the upper and lower bound for the 80% and 20% TSS of the range
    hold_ref = findBounds(ref)
    hold_tar = findBounds(tar)
    ref_upper = hold_ref[0]
    ref_lower = hold_ref[1]
    tar_upper = hold_tar[0]
    tar_lower = hold_tar[1]
    # Filters the reference dataset to 1) TSS >= 80% of the range and <= 20% of the range
    ref_eighty = filter_data(ref, '>=', ref_upper)
    ref_twenty = filter_data(ref, '<=', ref_lower)
    # Filters the target dataset to 1) TSS >= 80% of the range and <= 20% of the range
    tar_eighty = filter_data(tar, '>=', tar_upper)
    tar_twenty = filter_data(tar, '<=', tar_lower)
    
    
    
    # Matching the scores using the needleman wunsch algorithm with no alignments, minimum and maximum score
    # for each reference sequence is also gathered in this process.
    min_max_eighty = {}
    match_table_eighty = match(ref_eighty, tar_eighty, min_max_eighty)
    min_max_zero = {}
    match_table_zero = match(ref_twenty, tar_twenty, min_max_zero)
    
    # Calculate the percentage based on the match score and the min / max value of that reference sequence
    # Filters percentage matching >= 80% and <= 20%
    hold1 = filter_percentage(match_table_eighty, min_max_eighty, '>=', 0.8)
    hold2 = filter_percentage(match_table_eighty, min_max_eighty, '<=', 0.2)
    
    final_percentage_score_eighty_above_ninty = final_percentage_score_eighty_above_ninty.append(hold1)
    # final_percentage_score_zero_above_eighty = final_percentage_score_zero_above_eighty.append(filter_percentage(match_table_zero, min_max_zero, '>=', 0.9))
    final_percentage_score_eighty_below_twenty = final_percentage_score_eighty_below_twenty.append(hold2)
    # final_percentage_score_zero_below_twenty = final_percentage_score_zero_below_twenty.append(filter_percentage(match_table_zero, min_max_zero, '<=', 0.2))
    plot_chart(hold1, 50, 10, target + '_above_90.png')
    plot_chart(hold2, 50, 10, target + '_below_20.png')
    
# This function prepares the similarity matrix to the form that can be passed into the pairwise2 function
# csv(matrix) -> dict    
def prep_matrix(matrix):
    dict = {}
    for i, row in matrix.iterrows():
        row_name = matrix.iat[i,0]
        for j, row in matrix.iterrows():
            column_name = matrix.columns.values[j + 1]
            dict[row_name, column_name] = matrix.iat[i, j + 1];
    return dict
    
# filters data based on operator and target values
def filter_data(df, operator_string, target_value):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        if get_truth(df[df.columns[0]][i], operator_string, target_value):
            result= result.append(pd.DataFrame({df.columns.values[0]: df[df.columns[0]][i], 
                                             df.columns.values[1]: df[df.columns[1]][i]}, index=[0]), ignore_index=True)
    return result

# Needleman Wunsch Algorithm
# Also returns the min and max TSS score
def match(ref, tar, table):
    result = pd.DataFrame([])
    for i, row in ref.iterrows():
        print(i + 1, " / ", len(ref), " Done")
        ref_val = ref[ref.columns.values[1]][i]
        min_val = pairwise2.align.globalds(ref_val, tar[tar.columns.values[1]][0], matrix, -100, -100)[0][2]
        max_val = min_val
        for j, row in tar.iterrows():
            tar_val = tar[tar.columns.values[1]][j]
            score = pairwise2.align.globalds(ref_val, tar_val, matrix, -100, -100)[0][2]
            result = result.append(pd.DataFrame({'ref': ref_val, 
                                             'tar': tar_val, 'score': score}, index=[0]), ignore_index=True)
            min_val = min(min_val, score)
            max_val = max(max_val, score)
            table[ref_val] = (min_val, max_val)
    return result;

# Filter percentage based on operator and given percentage
def filter_percentage(df, bound, operator_string, target_percentage):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        ref_val = df[df.columns.values[0]][i]
        tar_val = df[df.columns.values[1]][i]
        min_score = bound[ref_val][0]
        max_score = bound[ref_val][1]
        total = max_score - min_score
        score = df[df.columns.values[2]][i]
        percentage = (score - min_score) / total
        if get_truth(percentage, operator_string, target_percentage):
            result = result.append(pd.DataFrame({'ref': ref_val, 
                                             'tar': tar_val, 'score': score, 
                                             'percentage': percentage, 'min_score_reference': min_score, 
                                             'max_score_reference': max_score}, index=[0]), 
        ignore_index=True)
    return result

# Given a dataset, returns the number that divides the upper 20% and the lower 80% of the range of data
# Given a dataset, returns the number that divides the lower 20% and the higher 80% of the range of data
def findBounds(df):
    min_TSS = min(df[df.columns.values[0]])
    max_TSS = max(df[df.columns.values[0]])
    lower_bound = min_TSS + 0.2 * (abs(min_TSS) + abs(max_TSS))
    upper_bound = min_TSS + 0.8 * (abs(min_TSS) + abs(max_TSS))
    return [upper_bound, lower_bound]

## Output plot
def plot_chart(df, num_1, num_2, name):
    fig = plt.figure(figsize=(num_1, num_2))
    plt.bar(df['ref'].value_counts().index, 
             df['ref'].value_counts())
    plt.xticks(rotation=70)
    plt.savefig(name ,bbox_inches='tight')
    plt.clf()

if __name__ == '__main__':
    compare_data_cross_species("Homo sapiens (Human)_Ameloblastin12_length.csv", "temp_f", "twelveAAHAmat.csv")
    