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

from contextlib import contextmanager

import biotite.structure as struc
import numpy as np


def report_stats(data: dict):
    """
    Generate and print a formatted table showing the statistics for each namespace in the given data.

    Args:
        data (dict): A dictionary where keys are namespaces and values are dictionaries of statistics.
    """
    columns = [list(data[namespace].keys()) for namespace in data]
    columns = sum(columns, [])
    columns = sorted(set(columns))
    output_lines = []
    output_lines.append(
        "                                     "
        + "".join([f"{c:>10s}" for c in columns])
    )
    for namespace in data:
        output_lines.append(
            f"{namespace[:35]:35s}  "
            + "".join([f"  {data[namespace].get(c, 0):8d}" for c in columns])
        )
    print("\n".join(output_lines))


class MMCIFStatsSummary:
    """
    A class for managing and summarizing MMCIF data statistics.
    """

    def __init__(self):
        self.logger_dict = {}
        self.data = {}

    def add_logger(self, pdbid: str) -> "MMCIFStatsLogger":
        """
        Add a new MMCIFStatsLogger instance to the logger dictionary for the given PDB ID.

        If a logger for the provided PDB ID already exists, a warning message will be printed,
        and a new logger will be initialized, overwriting the existing one.

        Args:
            pdbid (str): The PDB ID for which to add a logger.

        Returns:
            MMCIFStatsLogger: The newly created or overwritten logger instance for the given PDB ID.
        """
        if pdbid in self.logger_dict:
            print(f"WARNING: {pdbid} data stats logger exist, new logger initiated.")
        self.logger_dict[pdbid] = MMCIFStatsLogger(pdbid)
        return self.logger_dict[pdbid]

    def merge_data(self, data: dict):
        """
        Merge the provided data into the current instance's data.

        This method iterates over each namespace in the provided data and adds the
        key-value pairs to the corresponding namespace in the instance's data.
        If a namespace or key does not exist in the instance's data, it is created.

        Args:
            data (dict): A dictionary containing data to be merged.
                         The dictionary should have namespaces as keys and
                         another dictionary of key-value pairs as values.
        """
        for namespace in data:
            if namespace not in self.data:
                self.data[namespace] = {}
            for key in data[namespace]:
                if key not in self.data[namespace]:
                    self.data[namespace][key] = 0
                self.data[namespace][key] += int(data[namespace][key] != 0)

    def report(self):
        """
        Print a summary report of the statistics for all PDB IDs in the logger dictionary.

        This method generates a formatted table showing the statistics for each namespace
        across all PDB IDs, and prints it to the console.
        """
        print("----------------------------------------------------")
        print(f"statistics for all {len(self.logger_dict)} pdbids")
        report_stats(self.data)
        print("----------------------------------------------------")

    def get_data(self) -> dict:
        """
        Retrieve the data stored in the logger.

        This method returns the internal data dictionary that holds all the
        statistics and counts collected by the logger.

        Returns:
            dict: A dictionary containing the logged data.
        """
        return self.data

    def get_logger_data(self) -> dict:
        """
        Retrieve the data from all loggers in the logger dictionary.

        This method iterates over each logger in the `logger_dict` and calls the `get_data` method
        of each logger to obtain its data. It then creates a dictionary where the keys are the PDB IDs
        and the values are the corresponding logger data.

        Returns:
            dict: A dictionary containing the data from all loggers, with PDB IDs as keys.
        """
        return {pdbid: logger.get_data() for pdbid, logger in self.logger_dict.items()}

    def get_logger_dict(self) -> dict:
        """
        Retrieve the dictionary of loggers stored in the current instance.

        This method returns the `logger_dict` attribute, which is a dictionary
        containing MMCIFStatsLogger instances, with PDB IDs as keys.

        Returns:
            dict: A dictionary where keys are PDB IDs and values are MMCIFStatsLogger instances.
        """
        return self.logger_dict

    @classmethod
    def from_logger_dict(cls, logger_dict: dict) -> "MMCIFStatsSummary":
        """
        Create a new MMCIFStatsSummary instance from a dictionary of loggers.

        This class method initializes a new MMCIFStatsSummary object and populates it with
        loggers from the provided dictionary. It then merges the data from each logger into
        the summary's data attribute.

        Args:
            logger_dict (dict): A dictionary where keys are PDB IDs and values are MMCIFStatsLogger instances.

        Returns:
            MMCIFStatsSummary: A new MMCIFStatsSummary instance populated with the given loggers and their data.
        """
        summary = MMCIFStatsSummary()
        summary.logger_dict = logger_dict
        for logger in logger_dict.values():
            summary.merge_data(logger.data)
        return summary


class MMCIFStatsLogger:
    """
    A class for logging and tracking statistics related to MMCIF data.

    Args:
        pdbid (str): The PDB ID associated with the MMCIF data.
    """

    def __init__(self, pdbid: str = None):
        self.pdbid = pdbid
        self.data = {}
        self.atom_array = None

    def reset_pdbid(self, pdbid: str):
        """
        Reset the PDB ID associated with the MMCIFStatsLogger instance.

        This method updates the `pdbid` attribute of the logger with the provided value.

        Args:
            pdbid (str): The new PDB ID to be associated with the logger.
        """
        self.pdbid = pdbid

    def get_data(self) -> dict:
        """
        Retrieve the data stored in the logger.

        This method returns the internal data dictionary that holds all the
        statistics and counts collected by the logger.

        Returns:
            dict: A dictionary containing the logged data.
        """
        return self.data

    @contextmanager
    def log(self, namespace: str, atom_array: struc.AtomArray, components: list[str]):
        """
        Log the changes in component counts for a given namespace and atom array.

        This method is a context manager that records the initial counts of specified components
        in the AtomArray, yields the current instance, and then calculates and logs the changes
        in component counts after the context is exited.

        Args:
            namespace (str): The namespace under which the component counts will be logged.
            atom_array (AtomArray): The AtomArray containing the structural data.
            components (list[str]): A list of components to count, such as 'atom', 'bond', 'residue', etc.

        Yields:
            MMCIFStatsLogger: The current instance of the logger.
        """
        if namespace not in self.data:
            self.data[namespace] = {}
        self.atom_array = atom_array
        counts = self.count_component(components)
        yield self
        new_counts = self.count_component(components)
        for component, count, new_count in zip(components, counts, new_counts):
            if component not in self.data[namespace]:
                self.data[namespace][component] = 0
            self.data[namespace][component] += new_count - count
        self.atom_array = None

    def log_one(self, namespace: str, component: str, count: int = 1):
        """
        Log a single component count for a given namespace.

        This method is used to increment the count of a specific component within a given namespace.
        If the namespace or component does not exist in the data dictionary, it will be initialized.

        Args:
            namespace (str): The namespace under which the component count will be logged.
            component (str): The component to be logged.
            count (int, optional): The count to be added to the component. Defaults to 1.
        """
        if namespace not in self.data:
            self.data[namespace] = {}
        if component not in self.data[namespace]:
            self.data[namespace][component] = 0
        self.data[namespace][component] += int(count)

    def report(self):
        """
        Print a summary report of the statistics for the current PDB ID.

        This method generates a formatted table showing the statistics for each namespace
        in the current logger's data, and prints it to the console along with the PDB ID.
        """
        print("----------------------------------------------------")
        print(f"statistics for pdbid:{self.pdbid}")
        report_stats(self.data)
        print("----------------------------------------------------")

    def count_component(self, components: list[str] = None) -> list[int]:
        """
        Count the number of specified components in the current atom array.

        This method iterates through the provided list of components and counts
        the number of each component type in the atom array. If no components are
        provided, an empty list is used.

        Args:
            components (list[str], optional): A list of component types to count.
                Valid component types include 'atom', 'bond', 'residue', 'chain',
                'mol' and entity. Defaults to None.

        Returns:
            list[int]: A list of integers representing the counts of each component type.
        """
        assert self.atom_array is not None
        components = [] if components is None else components
        if len(self.atom_array) == 0:
            return [0] * len(components)
        counts = []
        for component_type in components:
            if component_type == "atom":
                counts.append(len(self.atom_array))
            elif component_type == "bond":
                counts.append(len(self.atom_array.bonds._bonds))
            elif component_type == "residue":
                counts.append(struc.get_residue_count(self.atom_array))
            elif component_type == "chain":
                counts.append(struc.get_chain_count(self.atom_array))
            elif component_type == "mol":
                counts.append(len(np.unique(self.atom_array.mol_id)))
            elif component_type == "entity":
                counts.append(len(np.unique(self.atom_array.label_entity_id)))
            else:
                raise NotImplementedError(f"Unknown component type: {component_type}")
        counts = [int(c) for c in counts]
        return counts
