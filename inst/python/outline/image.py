#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# =============================================================================
#      FileName: heline.py
#        Author: Chu Yanshuo
#         Email: yanshuochu@qq.com
# =============================================================================
"""


import numpy as np
from lib.base import Lines


class BorderLines(Lines):
    def __init__(self, pseudo_verts, verts):
        self._pseudo_verts = pseudo_verts
        self._verts = verts
        self._get_nbs_template()
        self._fake_spots = {}

    def _get_nbs_template(self):
        center = [2, 2]
        north = [2, 4]
        northeast = [3, 3]
        east = [4, 2]
        southeast = [3, 1]
        south = [2, 0]
        southwest = [1, 1]
        west = [0, 2]
        northwest = [1, 3]

        nbs_position = [
            north,
            northeast,
            east,
            southeast,
            south,
            southwest,
            west,
            northwest,
        ]
        nbs_index = [
            self._get_spot_index(self._pseudo_verts, nbp) for nbp in nbs_position
        ]
        nbs = [self._verts[index] for index in nbs_index]
        center_index = self._get_spot_index(self._pseudo_verts, center)
        offset = self._verts[center_index]
        self._nbs_template = np.array(nbs) - np.array(offset)

    def _get_neighbor_spots(self, pseudo_target_vert):
        pseudo_nb_d0 = [pseudo_target_vert[0], pseudo_target_vert[1] + 2]
        pseudo_nb_d1 = [
            pseudo_target_vert[0] + 1,
            pseudo_target_vert[1] + 1,
        ]
        pseudo_nb_d2 = [pseudo_target_vert[0] + 2, pseudo_target_vert[1]]
        pseudo_nb_d3 = [
            pseudo_target_vert[0] + 1,
            pseudo_target_vert[1] - 1,
        ]
        pseudo_nb_d4 = [pseudo_target_vert[0], pseudo_target_vert[1] - 2]
        pseudo_nb_d5 = [
            pseudo_target_vert[0] - 1,
            pseudo_target_vert[1] - 1,
        ]
        pseudo_nb_d6 = [pseudo_target_vert[0] - 2, pseudo_target_vert[1]]
        pseudo_nb_d7 = [
            pseudo_target_vert[0] - 1,
            pseudo_target_vert[1] + 1,
        ]

        pseudo_nbs = np.array(
            [
                pseudo_nb_d0,
                pseudo_nb_d1,
                pseudo_nb_d2,
                pseudo_nb_d3,
                pseudo_nb_d4,
                pseudo_nb_d5,
                pseudo_nb_d6,
                pseudo_nb_d7,
            ]
        )

        return self._map_neighbor_spots(pseudo_target_vert, pseudo_nbs)

    def _map_neighbor_spots(self, pseudo_target_vert, pseudo_nbs):
        nb_indices = [
            self._get_spot_index(self._pseudo_verts, point) for point in pseudo_nbs
        ]
        target_index = self._get_spot_index(self._pseudo_verts, pseudo_target_vert)
        target_vert = self._verts[target_index]
        nb_verts = []

        for d, index in enumerate(nb_indices):
            if index == -1:
                p_key = "{0}_{1}".format(pseudo_nbs[d][0], pseudo_nbs[d][1])
                if p_key in self._fake_spots:
                    temp_vert = self._fake_spots[p_key]
                else:
                    temp_vert = [
                        target_vert[0] + self._nbs_template[d][0],
                        target_vert[1] + self._nbs_template[d][1],
                    ]
                    self._fake_spots[p_key] = temp_vert
            else:
                temp_vert = self._verts[index]
            nb_verts.append(temp_vert)
        return nb_verts

    def _get_spot_index(self, verts, point):
        locs = (verts[:, 0] == point[0]) & (verts[:, 1] == point[1])
        if np.all(locs == False):
            target_index = -1
        else:
            target_index = np.where(locs)[0][0]

        return target_index

    def _get_hex_line_points(self, pseudo_target_vert, beta):
        target_index = self._get_spot_index(self._pseudo_verts, pseudo_target_vert)
        target_vert = self._verts[target_index]

        nbs = self._get_neighbor_spots(pseudo_target_vert)
        cross_points = self._get_cross_points(nbs, target_vert)

        hex_point_x_d0 = (cross_points[0][0] - target_vert[0]) * (
            2.0 / 3.0
        ) + target_vert[0]
        hex_point_y_d0 = (cross_points[0][1] - target_vert[1]) * (
            2.0 / 3.0
        ) + target_vert[1]

        hex_point_x_d1 = (nbs[1][0] - cross_points[1][0]) * (1.0 / 3.0) + cross_points[
            1
        ][0]
        hex_point_y_d1 = (nbs[1][1] - cross_points[1][1]) * (1.0 / 3.0) + cross_points[
            1
        ][1]

        hex_point_x_d3 = (nbs[3][0] - cross_points[1][0]) * (1.0 / 3.0) + cross_points[
            1
        ][0]
        hex_point_y_d3 = (nbs[3][1] - cross_points[1][1]) * (1.0 / 3.0) + cross_points[
            1
        ][1]

        hex_point_x_d4 = (cross_points[2][0] - target_vert[0]) * (
            2.0 / 3.0
        ) + target_vert[0]
        hex_point_y_d4 = (cross_points[2][1] - target_vert[1]) * (
            2.0 / 3.0
        ) + target_vert[1]

        hex_point_x_d5 = (nbs[5][0] - cross_points[3][0]) * (1.0 / 3.0) + cross_points[
            3
        ][0]
        hex_point_y_d5 = (nbs[5][1] - cross_points[3][1]) * (1.0 / 3.0) + cross_points[
            3
        ][1]

        hex_point_x_d7 = (nbs[7][0] - cross_points[3][0]) * (1.0 / 3.0) + cross_points[
            3
        ][0]
        hex_point_y_d7 = (nbs[7][1] - cross_points[3][1]) * (1.0 / 3.0) + cross_points[
            3
        ][1]

        hex_points = [
            [hex_point_x_d0, hex_point_y_d0],
            [hex_point_x_d1, hex_point_y_d1],
            [hex_point_x_d3, hex_point_y_d3],
            [hex_point_x_d4, hex_point_y_d4],
            [hex_point_x_d5, hex_point_y_d5],
            [hex_point_x_d7, hex_point_y_d7],
        ]

        return self._shrink_hex_points(hex_points, target_vert, beta)

    def _shrink_hex_points(self, nbs, target_vert, beta):
        return [self.__shrink_hex_point(nb, target_vert, beta) for nb in nbs]

    def __shrink_hex_point(self, point, target_vert, beta):
        nx = beta * (point[0] - target_vert[0]) + target_vert[0]
        ny = beta * (point[1] - target_vert[1]) + target_vert[1]
        return [nx, ny]

    def _get_cross_points(self, nbs, vert):
        cx_d0, cy_d0 = self.find_intersection(
            nbs[7][0],
            nbs[7][1],
            nbs[1][0],
            nbs[1][1],
            nbs[0][0],
            nbs[0][1],
            vert[0],
            vert[1],
        )
        cx_d2, cy_d2 = self.find_intersection(
            nbs[1][0],
            nbs[1][1],
            nbs[3][0],
            nbs[3][1],
            nbs[2][0],
            nbs[2][1],
            vert[0],
            vert[1],
        )
        cx_d4, cy_d4 = self.find_intersection(
            nbs[3][0],
            nbs[3][1],
            nbs[5][0],
            nbs[5][1],
            nbs[4][0],
            nbs[4][1],
            vert[0],
            vert[1],
        )
        cx_d6, cy_d6 = self.find_intersection(
            nbs[5][0],
            nbs[5][1],
            nbs[7][0],
            nbs[7][1],
            nbs[6][0],
            nbs[6][1],
            vert[0],
            vert[1],
        )

        return [[cx_d0, cy_d0], [cx_d2, cy_d2], [cx_d4, cy_d4], [cx_d6, cy_d6]]

    def get_borders(self, beta=1.0):
        points_array = [
            self._get_hex_line_points(ptv, beta) for ptv in self._pseudo_verts
        ]

        return [
            [
                [vert_nbs[0], vert_nbs[1]],
                [vert_nbs[1], vert_nbs[2]],
                [vert_nbs[2], vert_nbs[3]],
                [vert_nbs[4], vert_nbs[3]],
                [vert_nbs[5], vert_nbs[4]],
                [vert_nbs[0], vert_nbs[5]],
            ]
            for vert_nbs in points_array
        ]
