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

import matplotlib.pyplot as plt

import numpy as np

from itertools import islice

# DONE Only >= 80 now, implement <= 0
# DONE 1 --> Calculate upper and lower bound
# DONE 2 --> Get percentage, filter percentage (if statement) 
# Issue: Dict form of matrix
# 3 --> Plot

# This interface is used for the filter final percentage function
# final_percentage_score_zero_above_eighty['ref'].value_counts().plot('barh').invert_yaxis()

import operator

def get_truth(inp, relate, cut):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '=': operator.eq}
    return ops[relate](inp, cut)

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
    global final_percentage_score_eighty_above_eighty
    global final_percentage_score_zero_above_eighty
    global final_percentage_score_eighty_below_twenty
    global final_percentage_score_zero_below_twenty
    
    ref = pd.read_csv(reference)
    tar = pd.read_csv(target)
    
    matrix_unready = pd.read_csv(matrix_input)
    matrix = prep_matrix(matrix_unready)
    
    # Filtering
    ref_eighty = filter_data(ref, '>=', 80.0)
    ref_zero = filter_data(ref, '<=', 0.0)
    tar_eighty = filter_data(tar, '>=', 80.0)
    tar_zero = filter_data(tar, '<=', 0.0)
    
    # Matching
    min_max_eighty = {}
    match_table_eighty = match(ref_eighty, tar_eighty, min_max_eighty)
    min_max_zero = {}
    match_table_zero = match(ref_zero, tar_zero, min_max_zero)
    
    final_percentage_score_eighty_above_eighty = filter_percentage(match_table_eighty, min_max_eighty, '>=', 0.8)
    final_percentage_score_zero_above_eighty = filter_percentage(match_table_zero, min_max_zero, '>=', 0.8)
    final_percentage_score_eighty_below_twenty = filter_percentage(match_table_eighty, min_max_eighty, '<=', 0.2)
    final_percentage_score_zero_below_twenty = filter_percentage(match_table_zero, min_max_zero, '<=', 0.2)
    print("DONE")
    
def prep_matrix(matrix):
    dict = {}
    for i, row in matrix.iterrows():
        row_name = matrix.iat[i,0]
        for j, row in matrix.iterrows():
            column_name = matrix.columns.values[j + 1]
            dict[row_name, column_name] = matrix.iat[i, j + 1];
    return dict
    
def filter_data(df, operator_string, target_value):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        if get_truth(df[df.columns[0]][i], operator_string, target_value):
            result= result.append(pd.DataFrame({df.columns.values[0]: df[df.columns[0]][i], 
                                             df.columns.values[1]: df[df.columns[1]][i]}, index=[0]), ignore_index=True)
    return result

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

