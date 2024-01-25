#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: feature_plot.py
#        Author: Chu Yanshuo
#         Email: yanshuochu@qq.com
# =============================================================================
"""


import argparse
import os

import pandas as pd
import scanpy as sc

from feature_plot import FeaturePlot


class ModuleScorePlot(FeaturePlot):
    """Draw module score plot"""

    def __init__(self, adata_path, histology_annotation_path, ha_color_table_path=None):
        super().__init__(adata_path, histology_annotation_path, ha_color_table_path)

    def plot(
        self,
        marker_list_table_path,
        out_folder_path,
        to_plot_outline=True,
        to_plot_outline_type="pa",
    ):
        # step1: load marker list table #######################################
        marker_list_table = pd.read_csv(
            marker_list_table_path, header=0, na_filter=False, sep="\t"
        )
        for marker_list_name in marker_list_table.columns.tolist():
            marker_list = marker_list_table[marker_list_name].tolist()
            sc.tl.score_genes(self._adata, marker_list, score_name=marker_list_name)

            self._hex_figure_info_table = self._adata.obs[
                [
                    "x_array",
                    "y_array",
                    "Pathology_annotations",
                    "refined_pred",
                ]
            ]
            self._hex_figure_info_table["C"] = self._adata.obs[marker_list_name]

            self._hex_figure_info_adjust()

            self._make_plot(
                out_folder_path, marker_list_name, False, to_plot_outline_type
            )
            if to_plot_outline:
                if to_plot_outline_type == "all":
                    self._make_plot(out_folder_path, marker_list_name, True, "pa")
                    self._make_plot(out_folder_path, marker_list_name, True, "Cluster")
                else:
                    self._make_plot(
                        out_folder_path, marker_list_name, True, to_plot_outline_type
                    )


def main():
    parser = argparse.ArgumentParser(description="plot gene module score")

    parser.add_argument("--adataPath", dest="adataPath", help="adata file path")
    parser.add_argument(
        "--histologyAnnotationPath",
        dest="histologyAnnotationPath",
        help="histology annotation table path",
    )
    parser.add_argument(
        "--markerListTablePath",
        dest="markerListTablePath",
        help="marker list table path",
    )
    parser.add_argument(
        "--haColorTablePath",
        dest="haColorTablePath",
        default=None,
        help="histology annotation color table path",
    )
    parser.add_argument("--outFolderPath", dest="outFolderPath", help="out folder path")
    parser.add_argument(
        "--toPlotOutline", action="store_true", help="to plot outline or not"
    )
    parser.add_argument(
        "--toPlotOutlineType",
        dest="toPlotOutlineType",
        default="ha",
        help="outline type",
    )

    args = parser.parse_args()

    print(args)

    msp = ModuleScorePlot(
        adata_path=args.adataPath,
        histology_annotation_path=args.histologyAnnotationPath,
        ha_color_table_path=args.haColorTablePath,
    )
    msp.plot(
        marker_list_table_path=args.markerListTablePath,
        out_folder_path=args.outFolderPath,
        to_plot_outline=args.toPlotOutline,
        to_plot_outline_type=args.toPlotOutlineType,
    )


if __name__ == "__main__":
    main()
