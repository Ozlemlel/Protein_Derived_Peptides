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

# This function reads the data from the input argument
# Creates 2 global df, >=80 and <= 0 by calling filter_data_eighty and filter_data_zero
# also takes in target df
def read_data(reference, target):
    global ref_eighty
    global ref_zero
    global tar_eighty
    global tar_zero
    ref = pd.read_csv(reference)
    tar = pd.read_csv(target)
    ref_eighty = filter_data_eighty(ref)
    ref_zero = filter_data_zero(ref)
    tar_eighty = filter_data_eighty(tar)
    tar_zero = filter_data_zero(tar)
    
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

def match(ref, tar):
    result = pd.DataFrame([])
    for i, row in ref.iterrows():
        ref_val = ref[ref.columns.values[1]][i]
        for j, row in tar.iterrows():
            tar_val = tar[tar.columns.values[1]][j]
            score = pairwise2.align.globalmx(ref_val, tar_val, 1, -1)[0][2] / 12.0
            