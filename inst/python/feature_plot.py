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
import math
import os
import pickle

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from outline.hexbin import DrawOutline


class FeaturePlot(DrawOutline):
    """Draw feature plot"""

    def __init__(self, adata_path, histology_annotation_path, ha_color_table_path=None):
        self._adata_path = adata_path
        self._histology_annotation_path = histology_annotation_path
        self._ha_color_table_path = ha_color_table_path
        self._load_adata()

    def plot(
        self,
        marker_table_path,
        out_folder_path,
        to_plot_outline=True,
        to_plot_outline_type="pa",
    ):
        marker_table = pd.read_csv(
            marker_table_path, header=0, na_filter=False, sep="\t"
        )
        marker_table.columns.values[0] = "Marker"

        for marker in marker_table["Marker"].unique():
            self.plot_marker(
                marker, out_folder_path, to_plot_outline, to_plot_outline_type
            )

    def plot_marker(
        self, marker, out_folder_path, to_plot_outline=True, to_plot_outline_type="pa"
    ):
        if marker in self._adata.var_names:
            self._hex_figure_info_table = self._adata.obs[
                [
                    "x_array",
                    "y_array",
                    "Pathology_annotations",
                    "refined_pred",
                ]
            ]
            self._hex_figure_info_table["C"] = (
                self._adata[self._adata.obs.index, marker].X.toarray().flatten()
            )
            self._hex_figure_info_adjust()

            self._make_plot(out_folder_path, marker, False, to_plot_outline_type)
            if to_plot_outline:
                if to_plot_outline_type == "all":
                    self._make_plot(out_folder_path, marker, True, "pa")
                    self._make_plot(out_folder_path, marker, True, "Cluster")
                else:
                    self._make_plot(out_folder_path, marker, True, to_plot_outline_type)

    def _load_adata(self):
        # step1 merge histology annotations
        self._adata = sc.read(self._adata_path)
        histology_annotation_table = pd.read_csv(
            self._histology_annotation_path, header=0, na_filter=False, sep=","
        )
        histology_annotation_table.columns = ["Barcode", "Pathology_annotations"]
        if "Pathology_annotations" in self._adata.obs:
            self._adata.obs.drop("Pathology_annotations", axis=1)
        self._adata.obs = pd.merge(
            self._adata.obs,
            histology_annotation_table,
            left_index=True,
            right_on="Barcode",
            how="left",
        )
        self._adata.obs.set_index("Barcode", inplace=True)

    def _hex_figure_info_adjust(self):
        # swap x and y due to the hex shape
        self._hex_figure_info_table.columns = [
            "y_array",
            "x_array",
            "pathology_annotation",
            "cluster",
            "C",
        ]
        self._hex_figure_info_table["y_array"] = (
            self._hex_figure_info_table["y_array"].max()
            - self._hex_figure_info_table["y_array"]
        )
        self._hex_figure_info_table["x_array"] = (
            self._hex_figure_info_table["x_array"].max()
            - self._hex_figure_info_table["x_array"]
        )
        self._hex_figure_info_table.replace("", "Unlabelled", inplace=True)
        self._hex_figure_info_table.rename_axis("spot", inplace=True)

        self._unique_PAs = sorted(
            self._hex_figure_info_table["pathology_annotation"].unique()
        )
        mapping_array = pd.Series(
            [i + 1 for i in range(len(self._unique_PAs))], index=self._unique_PAs
        )
        self._hex_figure_info_table["pa_index"] = self._hex_figure_info_table[
            "pathology_annotation"
        ].map(mapping_array)

        self._unique_clusters = sorted(self._hex_figure_info_table["cluster"].unique())
        mapping_array = pd.Series(
            [i + 1 for i in range(len(self._unique_clusters))],
            index=self._unique_clusters,
        )
        self._hex_figure_info_table["cluster_index"] = self._hex_figure_info_table[
            "cluster"
        ].map(mapping_array)

    def _load_colors(self):
        # load colors #########################################################
        if self._ha_color_table_path is None:
            color_path_pa = os.path.join(
                os.path.dirname(self._adata_path),
                "pathology_annotation_colors.pkl",
            )
            if os.path.exists(color_path_pa):
                with open(color_path_pa, "rb") as color_file_pa:
                    self._color_cmap_pa = pickle.load(color_file_pa)
            else:
                self._color_cmap_pa = self._random_cmap(len(self._unique_PAs))
                with open(color_path_pa, "wb") as color_file_pa:
                    pickle.dump(self._color_cmap_pa, color_file_pa)
        else:
            self._ps_color_table = pd.read_csv(self._ha_color_table_path, sep="\t")
            color_values = [
                self._ps_color_table.loc[
                    self._ps_color_table["Tissue"] == t, "Color"
                ].iloc[0]
                for t in self._unique_PAs
            ]
            self._color_cmap_pa = mp.colors.ListedColormap(color_values)

        color_path_cluster = os.path.join(
            os.path.dirname(self._adata_path), "cluster_colors.pkl"
        )
        if os.path.exists(color_path_cluster):
            with open(color_path_cluster, "rb") as color_file_cluster:
                self._color_cmap_cluster = pickle.load(color_file_cluster)
        else:
            self._color_cmap_cluster = self._random_cmap(
                color_numbers=len(
                    np.unique(self._hex_figure_info_table["cluster_index"])
                )
            )
            with open(color_path_cluster, "wb") as color_file_cluster:
                pickle.dump(self._color_cmap_cluster, color_file_cluster)

    def _make_plot(
        self,
        out_folder_path,
        figure_name_prefix,
        to_plot_outline=True,
        to_plot_outline_type="pa",
    ):
        x_hb = self._hex_figure_info_table["x_array"]
        y_hb = self._hex_figure_info_table["y_array"]
        C_hb = self._hex_figure_info_table["C"]

        hb = plt.hexbin(
            x_hb,
            y_hb,
            C=C_hb,
            gridsize=(
                math.ceil(
                    (
                        max(self._hex_figure_info_table["x_array"])
                        - min(self._hex_figure_info_table["x_array"])
                    )
                    / 2
                ),
                math.ceil(
                    (
                        max(self._hex_figure_info_table["y_array"])
                        - min(self._hex_figure_info_table["y_array"])
                    )
                    / 2
                ),
            ),
            edgecolors="face",
            linewidth=0.0,
            cmap="gnuplot2",
            extent=(
                min(self._hex_figure_info_table["x_array"]) - 0.5,
                max(self._hex_figure_info_table["x_array"]) + 0.5,
                min(self._hex_figure_info_table["y_array"]) - 0.5,
                max(self._hex_figure_info_table["y_array"]) + 0.5,
            ),
        )
        plt.axis("off")
        cb = plt.colorbar(hb, shrink=0.6)
        cb.set_label("log1pEXP")
        plt.title(figure_name_prefix)
        plt.tight_layout()

        if to_plot_outline:
            self._load_colors()
            verts = hb.get_offsets()
            unique_y = np.unique(np.array([item[1] for item in verts]))
            diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
            unique_x = np.unique(np.array([item[0] for item in verts]))
            diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
            _spot_array = np.array(
                [
                    [x, y]
                    for x, y in zip(
                        self._hex_figure_info_table["x_array"],
                        self._hex_figure_info_table["y_array"],
                    )
                ]
            )
            index = self._get_nearest_index(verts, _spot_array)

            if to_plot_outline_type == "cluster" or to_plot_outline_type == "Cluster":
                cluster_indexed = self._hex_figure_info_table.index[index]
                cluster_index_df = self._hex_figure_info_table.loc[
                    cluster_indexed, "cluster_index"
                ]
                hex_colors_cluster = [
                    mp.colors.rgb2hex(color)
                    for color in self._color_cmap_cluster(
                        [
                            i
                            for i in range(
                                len(
                                    np.unique(
                                        self._hex_figure_info_table["cluster_index"]
                                    )
                                )
                            )
                        ]
                    )
                ]
                for tci, cid in enumerate(np.unique(cluster_index_df)):
                    if cid != "":
                        tcolor = hex_colors_cluster[tci]
                        self._plot_outline(
                            verts=verts,
                            pa=cluster_index_df,
                            tpa=cid,
                            plt=plt,
                            diff_y=diff_y,
                            diff_x=diff_x,
                            lw=0.5,
                            beta=0.88,
                            color=tcolor,
                        )
            elif to_plot_outline_type == "ha" or to_plot_outline_type == "pa":
                pa_indexed = self._hex_figure_info_table.index[index]
                pa_index_df = self._hex_figure_info_table.loc[pa_indexed, "pa_index"]

                hex_colors_pa = [
                    mp.colors.rgb2hex(color)
                    for color in self._color_cmap_pa(
                        [
                            i
                            for i in range(
                                len(np.unique(self._hex_figure_info_table["pa_index"]))
                            )
                        ]
                    )
                ]

                for tpi, tpa in enumerate(np.unique(pa_index_df)):
                    if tpa != "":
                        tcolor = hex_colors_pa[tpi]
                        self._plot_outline(
                            verts=verts,
                            pa=pa_index_df,
                            tpa=tpa,
                            plt=plt,
                            diff_y=diff_y,
                            diff_x=diff_x,
                            lw=0.5,
                            beta=0.88,
                            color=tcolor,
                        )
        if to_plot_outline:
            figure_path = "{0}/{1}_{2}.pdf".format(
                out_folder_path, figure_name_prefix, to_plot_outline_type
            )
        else:
            figure_path = "{0}/{1}.pdf".format(out_folder_path, figure_name_prefix)
        plt.savefig(figure_path)
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="plot gene expression")

    parser.add_argument("--adataPath", dest="adataPath", help="adata file path")
    parser.add_argument(
        "--histologyAnnotationPath",
        dest="histologyAnnotationPath",
        help="histology annotation table path",
    )
    parser.add_argument(
        "--markerTablePath", dest="markerTablePath", help="marker table to plot"
    )
    parser.add_argument(
        "--marker", dest="marker", default=None, help="marker"
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

    fp = FeaturePlot(
        adata_path=args.adataPath,
        histology_annotation_path=args.histologyAnnotationPath,
        ha_color_table_path=args.haColorTablePath,
    )
    if args.marker is None:
        fp.plot(
            marker_table_path=args.markerTablePath,
            out_folder_path=args.outFolderPath,
            to_plot_outline=args.toPlotOutline,
            to_plot_outline_type=args.toPlotOutlineType,
        )
    else:
        fp.plot(
            marker=args.marker,
            out_folder_path=args.outFolderPath,
            to_plot_outline=args.toPlotOutline,
            to_plot_outline_type=args.toPlotOutlineType,
        )


if __name__ == "__main__":
    main()
