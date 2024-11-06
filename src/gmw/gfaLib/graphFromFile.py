"Modelizes a graph object"
from itertools import count

from gfaLib.abstractions import GFALine, Orientation, reverse
from gfaLib.gfaparser import GFAParser
from gfaLib.abstractGraph import AbstractGraph

class GraphFromFile(AbstractGraph):

    def __init__(
        self,
        gfa_file: str | None = None,
        with_sequence: bool = True,
        low_memory: bool = False,
        regexp: str = ".*",
    ) -> None:
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
        self.paths[name] = {
            "id": identifier,
            "name": name,
            "origin": origin,
            "start_offset": start,
            "stop_offset": end if end is not None else sum([self.segments[x]['length'] for x, _ in chain]),
            "path": chain,
            **metadata
        }
