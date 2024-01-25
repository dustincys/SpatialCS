#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: plot_ha_cluster.py
#        Author: Chu Yanshuo
#         Email: yanshuoc@gmail.com
# =============================================================================
"""

import os

import numpy as np
from lib.base import Lines
from scipy.spatial.distance import cdist


class DrawOutline(Lines):
    """Functions to draw outline"""

    def _plot_outline(
        self, verts, pa, tpa, plt, diff_y, diff_x, lw, color="white", beta=0.88
    ):
        nv = verts[pa == tpa]
        lines = np.concatenate(
            [self._hex_lines(a=diff_y * 2 / 3, i=diff_x, off=off) for off in nv]
        )
        directions = np.concatenate([self._get_directions(off=off) for off in nv])
        uls, c = np.unique(lines.round(2), axis=0, return_counts=True)
        reshaped_array1 = uls[c == 1][:, np.newaxis]
        uls1_index = np.where(np.all(lines.round(2) == reshaped_array1, axis=(2, 3)))[1]
        t_d = directions[uls1_index]
        t_d = t_d[np.lexsort((t_d[:, 1], t_d[:, 0]))]

        unique_rows, indices = np.unique(t_d[:, :2], axis=0, return_inverse=True)

        grouped_values = np.split(t_d[:, 2], np.cumsum(np.bincount(indices))[:-1])
        grouped_values = [np.sort(np.unique(a)) for a in grouped_values]
        lines_s = np.concatenate(
            [
                self._hex_lines_direction(
                    a=diff_y * 2 / 3,
                    i=diff_x,
                    beta=beta - 0.1,
                    off=unique_rows[i],
                    directions=grouped_values[i].astype(int),
                )
                for i in range(len(unique_rows))
            ]
        )
        for l in lines_s:
            plt.plot(
                *l.transpose(), "w-", lw=lw, scalex=False, scaley=False, color="white"
            )
        lines_s = np.concatenate(
            [
                self._hex_lines_direction(
                    a=diff_y * 2 / 3,
                    i=diff_x,
                    beta=beta,
                    off=unique_rows[i],
                    directions=grouped_values[i].astype(int),
                )
                for i in range(len(unique_rows))
            ]
        )
        for l in lines_s:
            plt.plot(
                *l.transpose(), "w-", lw=lw, scalex=False, scaley=False, color=color
            )

    def _hex_lines(self, a=None, i=None, off=[0, 0]):
        """regular hexagon segment lines as `(xy1,xy2)` in clockwise
        order with points in line sorted top to bottom
        for irregular hexagon pass both `a` (vertical) and `i` (horizontal)"""
        if a is None:
            a = 2 / np.sqrt(3) * i
        if i is None:
            i = np.sqrt(3) / 2 * a
        h = a / 2
        xy = np.array(
            [
                [[0, a], [i, h]],
                [[i, h], [i, -h]],
                [[i, -h], [0, -a]],
                [[-i, -h], [0, -a]],  # flipped
                [[-i, h], [-i, -h]],  # flipped
                [[0, a], [-i, h]],  # flipped
            ]
        )
        return xy + off

    def _hex_lines_direction(
        self, a=None, i=None, beta=0.88, off=[0, 0], directions=[1, 2, 5]
    ):
        """regular hexagon segment lines as `(xy1,xy2)` in clockwise
        order with points in line sorted top to bottom
        for irregular hexagon pass both `a` (vertical) and `i` (horizontal)"""
        if a is None:
            a = 2 / np.sqrt(3) * i
        if i is None:
            i = np.sqrt(3) / 2 * a
        h = a / 2
        xy = np.array(
            [
                [[0, a], [i, h]],
                [[i, h], [i, -h]],
                [[i, -h], [0, -a]],
                [[-i, -h], [0, -a]],  # flipped
                [[-i, h], [-i, -h]],  # flipped
                [[0, a], [-i, h]],  # flipped
            ]
        )
        xy_beta = xy * beta
        xy_extents = np.array(
            [self._get_extent_points(i, xy, xy_beta) for i in range(len(xy_beta))]
        )
        for i_di, direction_i_di in enumerate(directions):
            i_di_left = i_di - 1
            i_di_right = (i_di + 1) % len(directions)
            if (directions[i_di_left] + 1) % 6 != direction_i_di:
                if direction_i_di < 3:
                    xy_beta[direction_i_di][0] = xy_extents[direction_i_di][0]
                else:
                    xy_beta[direction_i_di][1] = xy_extents[direction_i_di][0]
            if (directions[i_di_right] - 1) % 6 != direction_i_di:
                if direction_i_di < 3:
                    xy_beta[direction_i_di][1] = xy_extents[direction_i_di][1]
                else:
                    xy_beta[direction_i_di][0] = xy_extents[direction_i_di][1]
        xy_beta = xy_beta[directions]
        return xy_beta + off

    def _get_extent_points(self, i, xy, xy_beta):
        x1 = xy_beta[i][0][0]
        y1 = xy_beta[i][0][1]
        x2 = xy_beta[i][1][0]
        y2 = xy_beta[i][1][1]
        if i == 0:
            j = 5
        else:
            j = i - 1
        x3 = xy[j][0][0]
        y3 = xy[j][0][1]
        x4 = xy[j][1][0]
        y4 = xy[j][1][1]
        if i == 5:
            k = 0
        else:
            k = i + 1
        x5 = xy[k][0][0]
        y5 = xy[k][0][1]
        x6 = xy[k][1][0]
        y6 = xy[k][1][1]
        x_left, y_left = self.find_intersection(x1, y1, x2, y2, x3, y3, x4, y4)
        x_right, y_right = self.find_intersection(x1, y1, x2, y2, x5, y5, x6, y6)
        return [[x_left, y_left], [x_right, y_right]]


    def _get_directions(self, off=[0, 0]):
        xy = np.array(
            [
                [off[0], off[1], 0],
                [off[0], off[1], 1],
                [off[0], off[1], 2],
                [off[0], off[1], 3],
                [off[0], off[1], 4],
                [off[0], off[1], 5],
            ]
        )
        return xy

    def _get_nearest_index(self, array1, array2):
        distances = cdist(array1, array2)  # Calculate pairwise distances
        nearest_indices_c2 = np.argmin(distances, axis=1)
        return nearest_indices_c2

