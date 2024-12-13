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
    d = coordinate[..., None, :, :] - b[..., None, :]  #  [..., N_frame, N_atom, 3]
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
            [..., N_atom, 3]
        frame_atom_index (torch.Tensor): indices of three atoms in each frame
            [..., N_frame, 3] or [N_frame, 3]
        dim (torch.Tensor): along which dimension to select the frame atoms
    Returns:
        torch.Tensor: the constructed frames
            [..., N_frame, 3[three atom], 3[three coordinate]]
    """
    if len(frame_atom_index.shape) == 2:
        # the navie case
        x1 = torch.index_select(
            coordinate, dim=dim, index=frame_atom_index[:, 0]
        )  # [..., N_frame, 3]
        x2 = torch.index_select(
            coordinate, dim=dim, index=frame_atom_index[:, 1]
        )  # [..., N_frame, 3]
        x3 = torch.index_select(
            coordinate, dim=dim, index=frame_atom_index[:, 2]
        )  # [..., N_frame, 3]
        return torch.stack([x1, x2, x3], dim=dim)
    else:
        assert (
            frame_atom_index.shape[:dim] == coordinate.shape[:dim]
        ), "batch size dims should match"

    x1 = batched_gather(
        data=coordinate,
        inds=frame_atom_index[..., 0],
        dim=dim,
        no_batch_dims=len(coordinate.shape[:dim]),
    )  # [..., N_frame, 3]
    x2 = batched_gather(
        data=coordinate,
        inds=frame_atom_index[..., 1],
        dim=dim,
        no_batch_dims=len(coordinate.shape[:dim]),
    )  # [..., N_frame, 3]
    x3 = batched_gather(
        data=coordinate,
        inds=frame_atom_index[..., 2],
        dim=dim,
        no_batch_dims=len(coordinate.shape[:dim]),
    )  # [..., N_frame, 3]
    return torch.stack([x1, x2, x3], dim=dim)
