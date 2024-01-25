#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: plot_ha_cluster.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
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

from outline.hexbin import DrawOutline


class HACPlot(DrawOutline):
    """Draw combined plot"""

    def __init__(self, cluster_table_path, color_table_path=None):
        self._cluster_table_path = cluster_table_path
        self._color_table_path = color_table_path

        self._cluster_table = pd.read_csv(
            self._cluster_table_path,
            header=0,
            na_filter=False,
            index_col=0,
            sep="\t",
        )
        # swap x and y due to the hex shape
        self._cluster_table.columns = [
            "y_array",
            "x_array",
            "pathology_annotation",
            "cluster",
        ]
        self._cluster_table["y_array"] = (
            self._cluster_table["y_array"].max() - self._cluster_table["y_array"]
        )
        self._cluster_table["x_array"] = (
            self._cluster_table["x_array"].max() - self._cluster_table["x_array"]
        )
        self._cluster_table.replace("", "Unlabelled", inplace=True)
        self._cluster_table.rename_axis("spot", inplace=True)
        self._whr = (
            self._cluster_table["x_array"].max() - self._cluster_table["x_array"].min()
        ) / (
            self._cluster_table["y_array"].max() - self._cluster_table["y_array"].min()
        )

        if self._color_table_path is not None:
            self._color_table = pd.read_csv(self._color_table_path, sep="\t")

        self._unique_PAs = sorted(self._cluster_table["pathology_annotation"].unique())
        mapping_array = pd.Series(
            [i + 1 for i in range(len(self._unique_PAs))], index=self._unique_PAs
        )
        self._cluster_table["pa_index"] = self._cluster_table[
            "pathology_annotation"
        ].map(mapping_array)

        self._unique_clusters = sorted(self._cluster_table["cluster"].unique())
        mapping_array = pd.Series(
            [i + 1 for i in range(len(self._unique_clusters))],
            index=self._unique_clusters,
        )
        self._cluster_table["cluster_index"] = self._cluster_table["cluster"].map(
            mapping_array
        )

        self._load_colors()
        self._load_spot()

    def _load_spot(self):
        self._spot_array = np.array(
            [
                [x, y]
                for x, y in zip(
                    self._cluster_table["x_array"], self._cluster_table["y_array"]
                )
            ]
        )

    def _load_colors(self):
        # load colors #########################################################
        if self._color_table_path is None:
            color_path_pa = os.path.join(
                os.path.dirname(self._cluster_table_path),
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
            color_values = [
                self._color_table.loc[self._color_table["Tissue"] == t, "Color"].iloc[0]
                for t in self._unique_PAs
            ]
            self._color_cmap_pa = mp.colors.ListedColormap(color_values)

        color_path_cluster = os.path.join(
            os.path.dirname(self._cluster_table_path), "cluster_colors.pkl"
        )
        if os.path.exists(color_path_cluster):
            with open(color_path_cluster, "rb") as color_file_cluster:
                self._color_cmap_cluster = pickle.load(color_file_cluster)
        else:
            self._color_cmap_cluster = self._random_cmap(
                color_numbers=len(np.unique(self._cluster_table["cluster_index"]))
            )
            with open(color_path_cluster, "wb") as color_file_cluster:
                pickle.dump(self._color_cmap_cluster, color_file_cluster)

    def _dim_plot_ha_cluster(
        self,
        ax,
        ax_colorbar,
        highlight_PAs_fill=None,
        highlight_PAs_border=None,
        to_plot_color_bar=False,
        to_plot_outline=True,
        to_plot_outline_type="cluster",
    ):
        if highlight_PAs_fill is not None:
            x_hb = self._cluster_table["x_array"]
            y_hb = self._cluster_table["y_array"]

            if "NOTSHOW" not in highlight_PAs_fill:
                highlight_PAs_fill = highlight_PAs_fill + ["NOTSHOW"]

            mapping_array = pd.Series(
                [i + 1 for i in range(len(highlight_PAs_fill))], index=highlight_PAs_fill
            )

            self._cluster_table["pathology_annotation_temp"] = self._cluster_table[
                "pathology_annotation"
            ]
            self._cluster_table["pathology_annotation_temp"][
                ~np.isin(self._cluster_table["pathology_annotation"], highlight_PAs_fill)
            ] = "NOTSHOW"

            C_hb = self._cluster_table["pathology_annotation_temp"].map(mapping_array)

            color_values = [
                self._color_table.loc[self._color_table["Tissue"] == t, "Color"].iloc[0]
                for t in highlight_PAs_fill
            ]
            if "#FFFFFF" not in color_values:
                color_values = color_values + ["#FFFFFF"]

            temp_color_cmap_pa = mp.colors.ListedColormap(color_values)
        else:
            x_hb = self._cluster_table["x_array"]
            y_hb = self._cluster_table["y_array"]
            C_hb = self._cluster_table["pa_index"]
            temp_color_cmap_pa = self._color_cmap_pa

        hb = ax.hexbin(
            x_hb,
            y_hb,
            C=C_hb,
            gridsize=(
                math.ceil(
                    (
                        max(self._cluster_table["x_array"])
                        - min(self._cluster_table["x_array"])
                    )
                    / 2
                ),
                math.ceil(
                    (
                        max(self._cluster_table["y_array"])
                        - min(self._cluster_table["y_array"])
                    )
                    / 2
                ),
            ),
            edgecolors="face",
            linewidth=0.0,
            cmap=temp_color_cmap_pa,
            extent=(
                min(self._cluster_table["x_array"]) - 0.5,
                max(self._cluster_table["x_array"]) + 0.5,
                min(self._cluster_table["y_array"]) - 0.5,
                max(self._cluster_table["y_array"]) + 0.5,
            ),
        )
        ax.axis("off")
        ax.set_frame_on(False)

        if to_plot_color_bar:
            cb = plt.colorbar(hb, cax=ax_colorbar)
            cb.set_ticks(
                [
                    (i + 1.5) * (len(self._unique_PAs) - 1) / len(self._unique_PAs)
                    for i in range(len(self._unique_PAs))
                ]
            )
            cb.set_ticklabels(self._unique_PAs)

        if to_plot_outline:
            verts = hb.get_offsets()
            unique_y = np.unique(np.array([item[1] for item in verts]))
            diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
            unique_x = np.unique(np.array([item[0] for item in verts]))
            diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
            index = self._get_nearest_index(verts, self._spot_array)

            if to_plot_outline_type == "cluster":
                cluster_indexed = self._cluster_table.index[index]
                cluster_index_df = self._cluster_table.loc[
                    cluster_indexed, "cluster_index"
                ]
                hex_colors_cluster = [
                    mp.colors.rgb2hex(color)
                    for color in self._color_cmap_cluster(
                        [
                            i
                            for i in range(
                                len(np.unique(self._cluster_table["cluster_index"]))
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
                            plt=ax,
                            diff_y=diff_y,
                            diff_x=diff_x,
                            lw=0.5,
                            beta=0.88,
                            color=tcolor,
                        )
            else:
                pa_indexed = self._cluster_table.index[index]
                index_df = pd.Index(pa_indexed)
                pa = self._cluster_table.reindex(index_df)[
                    "pathology_annotation"
                ].str.strip()
                for tpa in np.unique(pa):
                    if tpa != "":
                        if highlight_PAs_border is not None:
                            if tpa not in highlight_PAs_border:
                                continue
                        tcolor = self._color_table.loc[
                            self._color_table["Tissue"] == tpa, "Color"
                        ].values[0]
                        self._plot_outline(
                            verts=verts,
                            pa=pa,
                            tpa=tpa,
                            plt=plt,
                            diff_y=diff_y,
                            diff_x=diff_x,
                            lw=0.5,
                            beta=0.88,
                            color=tcolor,
                        )

    def _dim_plot_cluster_ha(
        self,
        ax,
        ax_colorbar,
        highlight_PAs_fill=None,
        highlight_PAs_border=None,
        to_plot_color_bar=False,
        to_plot_outline=True,
        to_plot_outline_type="ha",
    ):
        if highlight_PAs_fill is not None:
            x_hb = self._cluster_table["x_array"]
            y_hb = self._cluster_table["y_array"]

            self._cluster_table["cluster_index_temp"] = self._cluster_table[
                "cluster_index"
            ]

            cluster_index_remain = self._cluster_table["cluster_index"][
                np.isin(self._cluster_table["pathology_annotation"], highlight_PAs_fill)
            ]
            self._cluster_table["cluster_index_temp"][
                ~np.isin(self._cluster_table["pathology_annotation"], highlight_PAs_fill)
            ] = -1

            mapping_array = pd.Series(
                [i + 1 for i in range(1 + len(np.unique(cluster_index_remain)))],
                index=np.append(np.unique(cluster_index_remain), -1),
            )
            C_hb = self._cluster_table["cluster_index_temp"].map(mapping_array)

            color_values = [
                mp.colors.rgb2hex(color)
                for color in self._color_cmap_cluster(
                    [ci - 1 for ci in np.unique(cluster_index_remain)]
                )
            ] + ["#FFFFFF"]
            temp_color_cmap_cluster = mp.colors.ListedColormap(color_values)
        else:
            x_hb = self._cluster_table["x_array"]
            y_hb = self._cluster_table["y_array"]
            C_hb = self._cluster_table["cluster_index"]
            temp_color_cmap_cluster = self._color_cmap_cluster

        hb = ax.hexbin(
            x_hb,
            y_hb,
            C=C_hb,
            gridsize=(
                math.ceil(
                    (
                        max(self._cluster_table["x_array"])
                        - min(self._cluster_table["x_array"])
                    )
                    / 2
                ),
                math.ceil(
                    (
                        max(self._cluster_table["y_array"])
                        - min(self._cluster_table["y_array"])
                    )
                    / 2
                ),
            ),
            edgecolors="face",
            linewidth=0.0,
            cmap=temp_color_cmap_cluster,
            extent=(
                min(self._cluster_table["x_array"]) - 0.5,
                max(self._cluster_table["x_array"]) + 0.5,
                min(self._cluster_table["y_array"]) - 0.5,
                max(self._cluster_table["y_array"]) + 0.5,
            ),
        )
        ax.axis("off")
        ax.set_frame_on(False)

        if to_plot_color_bar:
            cb = plt.colorbar(hb, cax=ax_colorbar)
            cb.set_ticks(
                [
                    (i + 1.5)
                    * (len(self._unique_clusters) - 1)
                    / len(self._unique_clusters)
                    for i in range(len(self._unique_clusters))
                ]
            )
            cb.set_ticklabels(self._unique_clusters)

        if to_plot_outline:
            verts = hb.get_offsets()
            unique_y = np.unique(np.array([item[1] for item in verts]))
            diff_y = np.sort(np.diff(np.sort(unique_y)))[0]
            unique_x = np.unique(np.array([item[0] for item in verts]))
            diff_x = np.sort(np.diff(np.sort(unique_x)))[0]
            index = self._get_nearest_index(verts, self._spot_array)

            if to_plot_outline_type == "ha":
                pa_indexed = self._cluster_table.index[index]
                index_df = pd.Index(pa_indexed)
                pa = self._cluster_table.reindex(index_df)[
                    "pathology_annotation"
                ].str.strip()
                for tpa in np.unique(pa):
                    if tpa != "":
                        if highlight_PAs_border is not None:
                            if tpa not in highlight_PAs_border:
                                continue
                        tcolor = self._color_table.loc[
                            self._color_table["Tissue"] == tpa, "Color"
                        ].values[0]
                        self._plot_outline(
                            verts=verts,
                            pa=pa,
                            tpa=tpa,
                            plt=plt,
                            diff_y=diff_y,
                            diff_x=diff_x,
                            lw=0.5,
                            beta=0.88,
                            color=tcolor,
                        )
            else:
                cluster_indexed = self._cluster_table.index[index]
                cluster_index_df = self._cluster_table.loc[
                    cluster_indexed, "cluster_index"
                ]
                hex_colors_cluster = [
                    mp.colors.rgb2hex(color)
                    for color in self._color_cmap_cluster(
                        [
                            i
                            for i in range(
                                len(np.unique(self._cluster_table["cluster_index"]))
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
                            plt=ax,
                            diff_y=diff_y,
                            diff_x=diff_x,
                            lw=0.5,
                            beta=0.88,
                            color=tcolor,
                        )

    def combined_dim_plot(
        self,
        highlight_PAs_fill=None,
        highlight_PAs_border=None,
        out_folder_path=".",
        out_file_name="spatial_hac.pdf",
    ):
        if highlight_PAs_fill is None:
            fig = plt.figure(
                figsize=(15, 15.0 / (0.7 * self._whr)), constrained_layout=True
            )
        else:
            fig = plt.figure(
                figsize=(15, 15.0 * 3 / (1.4 * self._whr)), constrained_layout=True
            )

        width_ratios = [0.1, 10, 0.1]
        spec = fig.add_gridspec(ncols=3, nrows=1, width_ratios=width_ratios)

        if highlight_PAs_fill is None:
            spec_main = spec[1].subgridspec(2, 2)
        else:
            spec_main = spec[1].subgridspec(3, 2)

        ax_left = fig.add_subplot(spec[0])
        ax_right = fig.add_subplot(spec[2])

        main_1 = fig.add_subplot(spec_main[0, 0])
        self._dim_plot_ha_cluster(
            main_1, ax_left, to_plot_color_bar=False, to_plot_outline=False
        )

        main_2 = fig.add_subplot(spec_main[1, 0])
        self._dim_plot_ha_cluster(
            main_2, ax_left, to_plot_color_bar=True, to_plot_outline=True
        )
        ax_left.yaxis.set_ticks_position("left")

        if highlight_PAs_fill is not None:
            main_3 = fig.add_subplot(spec_main[2, 0])
            self._dim_plot_ha_cluster(
                main_3,
                ax_left,
                highlight_PAs_fill=highlight_PAs_fill,
                highlight_PAs_border=highlight_PAs_border,
                to_plot_color_bar=False,
                to_plot_outline=True,
                to_plot_outline_type="ha",
            )

        main_4 = fig.add_subplot(spec_main[0, 1])
        self._dim_plot_cluster_ha(
            main_4, ax_right, to_plot_color_bar=False, to_plot_outline=False
        )

        main_5 = fig.add_subplot(spec_main[1, 1])
        self._dim_plot_cluster_ha(
            main_5, ax_right, to_plot_color_bar=True, to_plot_outline=True
        )
        ax_right.yaxis.set_ticks_position("right")

        if highlight_PAs_fill is not None:
            main_6 = fig.add_subplot(spec_main[2, 1])
            self._dim_plot_cluster_ha(
                main_6,
                ax_right,
                highlight_PAs_fill=highlight_PAs_fill,
                highlight_PAs_border=highlight_PAs_border,
                to_plot_color_bar=False,
                to_plot_outline=True,
                to_plot_outline_type="ha",
            )

        if not os.path.exists(out_folder_path):
            out_folder_path = "."
        out_file_path = os.path.join(out_folder_path, out_file_name)
        plt.savefig(out_file_path)
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="plot ha cluster")

    parser.add_argument(
        "--clusterTablePath", dest="clusterTablePath", help="cluster table path"
    )
    parser.add_argument(
        "--colorTablePath", dest="colorTablePath", help="color table path"
    )
    parser.add_argument("--outFolderPath", dest="outFolderPath", help="out folder path")
    parser.add_argument(
        "--outFileName",
        dest="outFileName",
        default="spatial_hac.pdf",
        help="out file name",
    )
    parser.add_argument(
        "--highlightFill",
        dest="highlightFill",
        nargs="+",
        default=None,
        help="clusters to show in fill color",
    )
    parser.add_argument(
        "--highlightBorder",
        dest="highlightBorder",
        nargs="+",
        default=None,
        help="clusters to show in border color",
    )
    args = parser.parse_args()

    print(args)

    if len(args.highlightFill) == 1 and args.highlightFill[0] == "":
        args.highlightFill = None

    if len(args.highlightBorder) == 1 and args.highlightBorder[0] == "":
        args.highlightBorder = None


    print(args)

    hac = HACPlot(
        cluster_table_path=args.clusterTablePath, color_table_path=args.colorTablePath
    )
    hac.combined_dim_plot(
        out_folder_path=args.outFolderPath,
        out_file_name=args.outFileName,
        highlight_PAs_fill=args.highlightFill,
        highlight_PAs_border=args.highlightBorder,
    )


if __name__ == "__main__":
    main()
