#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 22:05:05 2019

@author: Steve
"""

import pandas as pd

# This function reads the data from the input argument
# Creates 2 global df, >=80 and <= 0 by calling filter_data_eighty and filter_data_zero
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
        if df["TSSscoreStrong12"][i] >= 80.0:
            result= result.append(pd.DataFrame({'TSSscoreStrong12': df["TSSscoreStrong12"][i], 
                                             'twelvepeps': df["twelvepeps"][i]}, index=[0]), ignore_index=True)
    return result

def filter_data_zero(df):
    result = pd.DataFrame([])
    for i, row in df.iterrows():
        if df["TSSscoreStrong12"][i] <= 0.0:
            result = result.append(pd.DataFrame({'TSSscoreStrong12': df["TSSscoreStrong12"][i], 
                                             'twelvepeps': df["twelvepeps"][i]}, index=[0]), ignore_index=True)
    return result
