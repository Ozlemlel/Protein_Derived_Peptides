#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 22:05:05 2019

@author: Steve
"""

# Points?
# Percentage?

import pandas as pd

from Bio import pairwise2

from itertools import islice

# DONE Only >= 80 now, implement <= 0
# DONE 1 --> Calculate upper and lower bound
# 2 --> Get percentage, filter percentage (if statement) 
# 3 --> Plot

def read_data(reference, target):
    global ref_eighty
    global ref_zero
    global tar_eighty
    global tar_zero
    global match_table_eighty
    global match_table_zero
    global matrix
    global min_max_eighty
    global min_max_zero
    global final_percentage_score_eighty
    global final_percentage_score_zero
    
    matrix_unready = pd.read_csv('twelveAAHAmat.csv')
    matrix = prep_matrix(matrix_unready)
    
    ref = pd.read_csv(reference)
    tar = pd.read_csv(target)
    # Filtering
    ref_eighty = filter_data_eighty(ref)
    ref_zero = filter_data_zero(ref)
    tar_eighty = filter_data_eighty(tar)
    tar_zero = filter_data_zero(tar)
    # Matching
    min_max_eighty = {}
    match_table_eighty = match(ref_eighty, tar_eighty, min_max_eighty)
    min_max_zero = {}
    match_table_zero = match(ref_zero, tar_zero, min_max_zero)
    final_percentage_score_eighty = filter_percentage(match_table_eighty, min_max_eighty)
    final_percentage_score_zero = filter_percentage(match_table_zero, min_max_zero)
    
def prep_matrix(matrix):
    dict = {}
    for column_name, column in islice(matrix.transpose().iterrows(), 1, None):
        for i, row in matrix.iterrows():
            row_name = matrix[matrix.columns.values[0]][i]
            dict[row_name, column_name] = matrix[column_name][i];
    return dict
    
def filter_data_eighty(df):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        if df[df.columns[0]][i] >= 80.0:
            result= result.append(pd.DataFrame({df.columns.values[0]: df[df.columns[0]][i], 
                                             df.columns.values[1]: df[df.columns[1]][i]}, index=[0]), ignore_index=True)
    return result

def filter_data_zero(df):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        if df[df.columns[0]][i] <= 0.0:
            result = result.append(pd.DataFrame({df.columns.values[0]: df[df.columns[0]][i], 
                                             df.columns.values[1]: df[df.columns[1]][i]}, index=[0]), ignore_index=True)
    return result

def match(ref, tar, table):
    result = pd.DataFrame([])
    for i, row in ref.iterrows():
        print(i + 1, " / ", len(ref), " Done")
        ref_val = ref[ref.columns.values[1]][i]
        min_val = score = pairwise2.align.globaldx(ref_val, tar[tar.columns.values[1]][0], matrix)[0][2]
        max_val = min_val
        for j, row in tar.iterrows():
            tar_val = tar[tar.columns.values[1]][j]
            score = pairwise2.align.globaldx(ref_val, tar_val, matrix)[0][2]
            result = result.append(pd.DataFrame({'ref': ref_val, 
                                             'tar': tar_val, 'score': score}, index=[0]), ignore_index=True)
            min_val = min(min_val, score)
            max_val = max(max_val, score)
            table[ref_val] = (min_val, max_val)
    return result;

def filter_percentage(df, bound):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        ref_val = df[df.columns.values[0]][i]
        tar_val = df[df.columns.values[1]][i]
        min_score = bound[ref_val][0]
        max_score = bound[ref_val][1]
        total = max_score - min_score
        score = df[df.columns.values[2]][i]
        percentage = (score - min_score) / total
        if percentage >= 0.8:
            result = result.append(pd.DataFrame({'ref': ref_val, 
                                             'tar': tar_val, 'score': score, 'percentage': percentage}, index=[0]), ignore_index=True)
    return result

