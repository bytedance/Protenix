# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
from scipy.spatial.transform import Rotation


def angle_3p(a, b, c):
    """
    Calculate the angle between three points in a 2D space.

    Args:
        a (list or array-like): The coordinates of the first point.
        b (list or array-like): The coordinates of the second point.
        c (list or array-like): The coordinates of the third point.

    Returns:
        float: The angle in degrees (0, 180) between the vectors
               from point a to point b and point b to point c.
    """
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    ab = b - a
    bc = c - b

    dot_product = np.dot(ab, bc)

    norm_ab = np.linalg.norm(ab)
    norm_bc = np.linalg.norm(bc)

    cos_theta = np.clip(dot_product / (norm_ab * norm_bc + 1e-4), -1, 1)
    theta_radians = np.arccos(cos_theta)
    theta_degrees = np.degrees(theta_radians)
    return theta_degrees


def random_transform(
    points, max_translation=1.0, apply_augmentation=False, centralize=True
) -> np.ndarray:
    """
    Randomly transform a set of 3D points.

    Args:
        points (numpy.ndarray): The points to be transformed, shape=(N, 3)
        max_translation (float): The maximum translation value. Default is 1.0.
        apply_augmentation (bool): Whether to apply random rotation/translation on ref_pos

    Returns:
        numpy.ndarray: The transformed points.
    """
    if centralize:
        points = points - points.mean(axis=0)
    if not apply_augmentation:
        return points
    translation = np.random.uniform(-max_translation, max_translation, size=3)
    R = Rotation.random().as_matrix()
    transformed_points = np.dot(points + translation, R.T)
    return transformed_points
