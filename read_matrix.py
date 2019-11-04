#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 00:08:32 2019

@author: mabochen
"""

import pandas as pd
from itertools import islice

matrix = pd.read_csv('twelveAAHAmat.csv')

dict = {}

for column_name, column in islice(matrix.transpose().iterrows(), 1, None):
    for i, row in matrix.iterrows():
        row_name = matrix[matrix.columns.values[0]][i]
        dict[row_name, column_name] = matrix[column_name][i];