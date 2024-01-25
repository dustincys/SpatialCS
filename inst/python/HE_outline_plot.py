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
import json
import os

import matplotlib as mp
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib.base import load_colors
from outline.image import BorderLines


class HEOutlinePlot:
    """Draw module score plot"""

    def __init__(
        self,
        spatialranger_out_path,
        cluster_table_path,
        ha_color_table_path=None,
    ):
        self._spatialranger_out_path = spatialranger_out_path
        self._cluster_table_path = cluster_table_path
        self._ha_color_table_path = ha_color_table_path
        self._load_adata()
        self._load_spatial()
        self._load_image()

    def _load_adata(self):
        self._cluster_table = pd.read_csv(
            self._cluster_table_path,
            header=0,
            na_filter=False,
            index_col=0,
            sep="\t",
        )
        self._cluster_table.columns = [
            "x_array",
            "y_array",
            "pathology_annotation",
            "cluster",
        ]

    def _load_spatial(self):
        # step1: read json scale
        with open(
            os.path.join(
                self._spatialranger_out_path, "spatial", "scalefactors_json.json"
            ),
            "r",
        ) as json_file:
            scalefactors_dict = json.load(json_file)
            tissue_hires_scalef = float(scalefactors_dict["tissue_hires_scalef"])

        # step2: read tissue positions
        tissue_positions = pd.read_csv(
            os.path.join(
                self._spatialranger_out_path, "spatial", "tissue_positions.csv"
            )
        )

        tissue_positions["pxl_row_in_hires"] = (
            tissue_hires_scalef * tissue_positions["pxl_row_in_fullres"]
        )
        tissue_positions["pxl_col_in_hires"] = (
            tissue_hires_scalef * tissue_positions["pxl_col_in_fullres"]
        )
        tissue_positions = tissue_positions.set_index("barcode")
        self._all_spots_info = pd.merge(
            self._cluster_table,
            tissue_positions,
            left_index=True,
            right_index=True,
            how="right",
        )

    def _load_image(self):
        self._img = mpimg.imread(
            os.path.join(
                self._spatialranger_out_path, "spatial", "tissue_hires_image.png"
            )
        )

    def _get_borders(self, beta=0.8):
        pseudo_verts = np.array(
            [
                [x, y]
                for y, x in zip(
                    self._all_spots_info["array_row"],
                    self._all_spots_info["array_col"],
                )
            ]
        )
        verts = np.array(
            [
                [x, y]
                for y, x in zip(
                    self._all_spots_info["pxl_row_in_hires"],
                    self._all_spots_info["pxl_col_in_hires"],
                )
            ]
        )
        img_b = BorderLines(pseudo_verts, verts)

        return np.array(img_b.get_borders(beta))

    def plot(
        self,
        out_folder_path,
        figure_name_prefix,
        to_plot_outline_type="pa",
    ):
        all_borders = self._get_borders(beta = 1.0)
        all_borders_beta = self._get_borders(beta = 0.8)
        tissue_spots, tissue_spots_indices = self._get_tissue_spots()
        tissue_spots = tissue_spots.reset_index()
        tissue_spot_borders = all_borders[tissue_spots_indices]
        tissue_spot_borders_beta = all_borders_beta[tissue_spots_indices]

        colors_ha, colors_cluster = self._get_colors(tissue_spots)

        plt.imshow(self._img, origin="upper")

        if to_plot_outline_type == "ha" or to_plot_outline_type == "pa":
            for color_i, pai in enumerate(sorted(tissue_spots["pa_index"].unique())):
                pai_index = tissue_spots.index[
                    tissue_spots["pa_index"] == pai
                ]

                lines = np.concatenate(tissue_spot_borders[pai_index])
                lines_beta = np.concatenate(tissue_spot_borders_beta[pai_index])
                uls, c = np.unique(lines.round(2), axis=0, return_counts=True)
                reshaped_array1 = uls[c == 1][:, np.newaxis]
                uls1_index = np.where(
                    np.all(lines.round(2) == reshaped_array1, axis=(2, 3))
                )[1]
                lines_beta = lines_beta[uls1_index]
                lines = lines[uls1_index]

                for l_beta in lines_beta:
                    plt.plot(
                        *l_beta.transpose(), "w-", lw=0.5, scalex=False, scaley=False, color=colors_ha[color_i]
                    )

                for l in lines:
                    plt.plot(
                        *l.transpose(), "w-", lw=0.1, scalex=False, scaley=False, color="white"
                    )

        else:
            for color_i, cai in enumerate(sorted(tissue_spots["cluster_index"].unique())):
                cai_index = tissue_spots.index[
                    tissue_spots["cluster_index"] == cai
                ]

                lines = np.concatenate(tissue_spot_borders[cai_index])
                lines_beta = np.concatenate(tissue_spot_borders_beta[cai_index])
                uls, c = np.unique(lines.round(2), axis=0, return_counts=True)
                reshaped_array1 = uls[c == 1][:, np.newaxis]
                uls1_index = np.where(
                    np.all(lines.round(2) == reshaped_array1, axis=(2, 3))
                )[1]

                lines_beta = lines_beta[uls1_index]
                lines = lines[uls1_index]

                for l_beta in lines_beta:
                    plt.plot(
                        *l_beta.transpose(), "w-", lw=0.5, scalex=False, scaley=False, color=colors_cluster[color_i]
                    )
                for l in lines:
                    plt.plot(
                        *l.transpose(), "w-", lw=0.1, scalex=False, scaley=False, color="white"
                    )

        figure_path = "{0}/{1}_{2}.pdf".format(
            out_folder_path, figure_name_prefix, to_plot_outline_type
        )
        plt.axis("off")
        plt.savefig(figure_path)
        plt.close()

    def _get_colors(self, tissue_spots):
        unique_PAs = sorted(tissue_spots["pathology_annotation"].unique())
        cluster_num = len(sorted(tissue_spots["cluster"].unique()))
        color_cmap_ha, color_cmap_cluster = load_colors(
            self._cluster_table_path, self._ha_color_table_path, unique_PAs, cluster_num
        )

        hex_colors_ha = [
            mp.colors.rgb2hex(color)
            for color in color_cmap_ha(
                [i for i in range(len(np.unique(tissue_spots["pa_index"])))]
            )
        ]

        hex_colors_cluster = [
            mp.colors.rgb2hex(color)
            for color in color_cmap_cluster(
                [i for i in range(len(np.unique(tissue_spots["cluster_index"])))]
            )
        ]

        return hex_colors_ha, hex_colors_cluster

    def _get_tissue_spots(self):
        tissue_spots = self._all_spots_info[self._all_spots_info["in_tissue"] == 1]
        tissue_spots_indices = np.array(
            [
                index
                for index, value in enumerate(self._all_spots_info.index)
                if value in tissue_spots.index
            ]
        )

        tissue_spots = tissue_spots[
            [
                "array_row",
                "array_col",
                "pxl_row_in_hires",
                "pxl_col_in_hires",
                "pathology_annotation",
                "cluster",
            ]
        ]

        tissue_spots.replace("", "Unlabelled", inplace=True)
        tissue_spots.rename_axis("spot", inplace=True)

        unique_PAs = sorted(tissue_spots["pathology_annotation"].unique())
        mapping_array = pd.Series(
            [i + 1 for i in range(len(unique_PAs))], index=unique_PAs
        )
        tissue_spots["pa_index"] = tissue_spots["pathology_annotation"].map(
            mapping_array
        )

        unique_clusters = sorted(tissue_spots["cluster"].unique())
        mapping_array = pd.Series(
            [i + 1 for i in range(len(unique_clusters))],
            index=unique_clusters,
        )
        tissue_spots["cluster_index"] = tissue_spots["cluster"].map(mapping_array)

        return tissue_spots, tissue_spots_indices


def main():
    parser = argparse.ArgumentParser(description="plot HE outline")

    parser.add_argument(
        "--spatialrangerOutPath",
        dest="spatialrangerOutPath",
        help="spatial ranger output path",
    )
    parser.add_argument(
        "--clusterTablePath",
        dest="clusterTablePath",
        help="cluster table path",
    )
    parser.add_argument(
        "--histologyAnnotationPath",
        dest="histologyAnnotationPath",
        help="histology annotation table path",
    )
    parser.add_argument(
        "--haColorTablePath",
        dest="haColorTablePath",
        default=None,
        help="histology annotation color table path",
    )
    parser.add_argument("--outFolderPath", dest="outFolderPath", help="out folder path")
    parser.add_argument(
        "--figureNamePrefix",
        dest="figureNamePrefix",
        default="hexbin_outline",
        help="hexbin outline",
    )
    parser.add_argument(
        "--toPlotOutlineType",
        dest="toPlotOutlineType",
        default="ha",
        help="outline type",
    )

    args = parser.parse_args()

    print(args)

    heop = HEOutlinePlot(
        spatialranger_out_path=args.spatialrangerOutPath,
        cluster_table_path=args.clusterTablePath,
        ha_color_table_path=args.haColorTablePath,
    )

    heop.plot(
        out_folder_path=args.outFolderPath,
        figure_name_prefix=args.figureNamePrefix,
        to_plot_outline_type=args.toPlotOutlineType,
    )


if __name__ == "__main__":
    main()
