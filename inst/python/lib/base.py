#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: base.py
#        Author: Chu Yanshuo
#         Email: yanshuochu@qq.com
# =============================================================================
"""

import os
import pickle

import matplotlib as mp
import numpy as np
import pandas as pd
from distinctipy import distinctipy


class Lines:
    def find_intersection(self, x1, y1, x2, y2, x3, y3, x4, y4):
        if np.isinf((y2 - y1) / (x2 - x1)):
            # Line 1 is vertical
            x_intersect = x1
            slope2 = (y4 - y3) / (x4 - x3)
            intercept2 = y3 - slope2 * x3
            y_intersect = slope2 * x_intersect + intercept2
        elif np.isinf((y4 - y3) / (x4 - x3)):
            # Line 2 is vertical
            x_intersect = x3
            slope1 = (y2 - y1) / (x2 - x1)
            intercept1 = y1 - slope1 * x1
            y_intersect = slope1 * x_intersect + intercept1
        else:
            # Neither line is vertical, calculate intersection normally
            slope1 = (y2 - y1) / (x2 - x1)
            slope2 = (y4 - y3) / (x4 - x3)
            intercept1 = y1 - slope1 * x1
            intercept2 = y3 - slope2 * x3
            x_intersect = (intercept2 - intercept1) / (slope1 - slope2)
            y_intersect = slope1 * x_intersect + intercept1
        # Return the coordinates of the intersection point
        return x_intersect, y_intersect

def random_cmap(color_numbers=10):
    color_tuples = distinctipy.get_colors(color_numbers)
    return mp.colors.ListedColormap(color_tuples)

def load_colors(adata_path, ha_color_table_path, unique_HAs, cluster_num):
    color_cmap_pa, color_cmap_cluster = None, None
    if ha_color_table_path is None:
        color_path_pa = os.path.join(
            os.path.dirname(adata_path), "pathology_annotation_colors.pkl"
        )
        if os.path.exists(color_path_pa):
            with open(color_path_pa, "rb") as color_file_pa:
                color_cmap_pa = pickle.load(color_file_pa)
        else:
            color_cmap_pa = random_cmap(len(unique_HAs))
            with open(color_path_pa, "wb") as color_file_pa:
                pickle.dump(color_cmap_pa, color_file_pa)
    else:
        ps_color_table = pd.read_csv(ha_color_table_path, sep="\t")
        color_values = [
            ps_color_table.loc[ps_color_table["Tissue"] == t, "Color"].iloc[0]
            for t in unique_HAs
        ]
        color_cmap_pa = mp.colors.ListedColormap(color_values)

    color_path_cluster = os.path.join(os.path.dirname(adata_path), "cluster_colors.pkl")
    if os.path.exists(color_path_cluster):
        with open(color_path_cluster, "rb") as color_file_cluster:
            color_cmap_cluster = pickle.load(color_file_cluster)
    else:
        color_cmap_cluster = random_cmap(color_numbers=cluster_num)
        with open(color_path_cluster, "wb") as color_file_cluster:
            pickle.dump(color_cmap_cluster, color_file_cluster)

    return color_cmap_pa, color_cmap_cluster
