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

import torch
import torch.nn.functional as F

from protenix.model.utils import batched_gather


def expressCoordinatesInFrame(
    coordinate: torch.Tensor, frames: torch.Tensor, eps: float = 1e-8
) -> torch.Tensor:
    """Algorithm 29 Express coordinate in frame

    Args:
        coordinate (torch.Tensor): the input coordinate
            [..., N_atom, 3]
        frames (torch.Tensor): the input frames
            [..., N_frame, 3, 3]
        eps (float): Small epsilon value

    Returns:
        torch.Tensor: the transformed coordinate projected onto frame basis
            [..., N_frame, N_atom, 3]
    """
    # Extract frame atoms
    a, b, c = torch.unbind(frames, dim=-2)  # a, b, c shape: [..., N_frame, 3]
    w1 = F.normalize(a - b, dim=-1, eps=eps)
    w2 = F.normalize(c - b, dim=-1, eps=eps)
    # Build orthonormal basis
    e1 = F.normalize(w1 + w2, dim=-1, eps=eps)
    e2 = F.normalize(w2 - w1, dim=-1, eps=eps)
    e3 = torch.cross(e1, e2, dim=-1)  # [..., N_frame, 3]
    # Project onto frame basis
    d = coordinate[..., None, :, :] - b[..., None, :]  # [..., N_frame, N_atom, 3]
    x_transformed = torch.cat(
        [
            torch.sum(d * e1[..., None, :], dim=-1, keepdim=True),
            torch.sum(d * e2[..., None, :], dim=-1, keepdim=True),
            torch.sum(d * e3[..., None, :], dim=-1, keepdim=True),
        ],
        dim=-1,
    )  # [..., N_frame, N_atom, 3]
    return x_transformed


def gather_frame_atom_by_indices(
    coordinate: torch.Tensor, frame_atom_index: torch.Tensor, dim: int = -2
) -> torch.Tensor:
    """construct frames from coordinate

    Args:
        coordinate (torch.Tensor):  the input coordinate
            [..., N_atom, 3[three coordinates]]
        frame_atom_index (torch.Tensor): indices of three atoms in each frame
            [..., N_frame, 3[three atoms per frame]] or [N_frame, 3[three atoms per frame]]
        dim (int): along which dimension to select the frame atoms
    Returns:
        torch.Tensor: the constructed frames
            [..., N_frame, 3[three atoms per frame], 3[three coordinates]]
    """
    if len(frame_atom_index.shape) == 2:
        # the naive case
        return coordinate[..., frame_atom_index, :]
    else:
        assert (
            frame_atom_index.shape[:dim] == coordinate.shape[:dim]
        ), f"the size of each batch dim should match, got {frame_atom_index.shape[:dim]} and {coordinate.shape[:dim]}"

    reshaped_frame_atom_index = frame_atom_index.reshape(*frame_atom_index.shape[:-2], -1)  # [..., N_frame*3]
    batched_frame_atom_coordinates = batched_gather(
        data=coordinate,
        inds=reshaped_frame_atom_index
    )  # [..., N_frame*3, 3[three coordinates]]
    return batched_frame_atom_coordinates.reshape(*batched_frame_atom_coordinates.shape[:-2], frame_atom_index.shape[-2], frame_atom_index.shape[-1], coordinate.shape[-1])  # [..., N_frame, 3, 3]
