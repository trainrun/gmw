""" 
GMW: Genomic Microbe-Wise - hybrid assembly and contamination removal tool 

Copyright (C) 2025 Wenbing Chen 
www.github.com/trainrun/gmw 

License: 
This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   
See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see <https://www.gnu.org/licenses/>. 
"""

from itertools import count

from gfaLib.abstractions import GFALine, Orientation, reverse
from gfaLib.gfaparser import GFAParser
from gfaLib.abstractGraph import AbstractGraph

class GraphFromFile(AbstractGraph):
    """
    A class for creating a GFA graph from a GFA file.
    
    This class reads a GFA (Graphical Fragment Assembly) file and constructs
    a graph representation with segments, paths, links, and headers. It supports
    different memory modes and can handle various GFA formats.
    """

    def __init__(
        self,
        gfa_file: str | None = None,
        with_sequence: bool = True,
        low_memory: bool = False,
        regexp: str = ".*",
    ) -> None:
        """
        Initialize a GraphFromFile object by parsing a GFA file.
        
        Parameters:
            gfa_file: Path to the GFA file to parse
            with_sequence: Whether to load sequence data into memory
            low_memory: If True, uses memory-efficient mode with reduced functionality
            regexp: Regular expression pattern to filter lines during parsing
        """
        super().__init__()
        # Declaring format attributes, generators...
        self.metadata: dict = {
            'version': GFAParser.get_gfa_format(gfa_file_path=gfa_file) if gfa_file and not low_memory else 'unknown',
            'next_node_name': (x for x in count(start=1) if str(x) not in self.segments) if not low_memory else 'unknown',
            'with_sequence': with_sequence
        }

        with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
            for gfa_line in gfa_reader:

                name, line_type, datas = GFAParser.read_gfa_line(
                    datas=[__.strip() for __ in gfa_line.split('\t')],
                    load_sequence_in_memory=with_sequence and not low_memory,
                    regexp_pattern=regexp,
                    memory_mode=low_memory,
                )
                if line_type == GFALine.SEGMENT:
                    self.segments[name] = datas
                elif line_type in (GFALine.WALK, GFALine.PATH):
                    self.paths[name] = datas
                elif line_type == GFALine.LINK:
                    while name in self.lines:
                        name = (name[0], name[1], name[2] + 1)
                    self.lines[name] = datas
                    # if name not in self.lines:
                    #     self.lines[name] = datas
                    # else:
                    #     [_ors,] = datas["orientation"]
                    #     self.lines[name]["orientation"].add((_ors[0], _ors[1]))

                elif line_type == GFALine.HEADER:
                    self.headers.append(datas)
                else:
                    pass  # Ignore unknown line types

    def add_node(self, name, sequence, **metadata):
        """
        Add a new node (segment) to the graph.
        
        Parameters:
            name: Unique identifier for the node
            sequence: The DNA/RNA sequence of the segment
            **metadata: Additional attributes to store with the segment
        """
        if not self.metadata['with_sequence']:
            self.segments[name] = {
                'length': len(sequence),
                **metadata
            }
        else:
            self.segments[name] = {
                'seq': sequence,
                'length': len(sequence),
                **metadata
            }

    def add_edge(self, source, ori_source, sink, ori_sink, **metadata):
        """
        Add a new edge (link) between two nodes in the graph.
        
        Parameters:
            source: Name of the source node
            ori_source: Orientation of the source node ('+', '-', '?', '=')
            sink: Name of the sink node
            ori_sink: Orientation of the sink node ('+', '-', '?', '=')
            **metadata: Additional attributes to store with the edge
            
        Raises:
            ValueError: If the orientation values are not compatible with GFA format
        """
        if not ori_sink in ['+', '-', '?', '=']:
            try:
                ori_sink = ori_sink.value
            except:
                raise ValueError("Not compatible with GFA format.")
        if not ori_source in ['+', '-', '?', '=']:
            try:
                ori_source = ori_source.value
            except:
                raise ValueError("Not compatible with GFA format.")
        if (source, sink) not in self.lines:
            self.lines[(source, sink)] = {
                'orientation': set([(Orientation(ori_source), Orientation(ori_sink))]),
                **metadata
            }
        else:
            self.lines[(source, sink)]['orientation'] = self.lines[(source, sink)].get(
                'orientation', set()) | set([(Orientation(ori_source), Orientation(ori_sink))])

    def add_path(self, identifier, name, chain, start=0, end=None, origin=None, **metadata):
        """
        Add a path (ordered sequence of nodes) to the graph.
        
        Parameters:
            identifier: Unique identifier for the path
            name: Name of the path
            chain: List of (node_name, orientation) tuples defining the path
            start: Starting offset in the path (default: 0)
            end: Ending offset in the path (default: sum of all segment lengths)
            origin: Source or origin information for the path
            **metadata: Additional attributes to store with the path
        """
        self.paths[name] = {
            "id": identifier,
            "name": name,
            "origin": origin,
            "start_offset": start,
            "stop_offset": end if end is not None else sum([self.segments[x]['length'] for x, _ in chain]),
            "path": chain,
            **metadata
        }
