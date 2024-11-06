from gfaLib.abstractions import Orientation
from gfaLib.gfaparser import GFAParser

class AbstractGraph:

    def __init__(self):
        self.segments: dict[str, dict] = {}
        self.lines: dict[tuple[str, str, int], dict] = {}
        self.paths: dict[str, dict] = {}
        self.headers: list[dict] = []
        
    def save_graph(self, output_file, minimal=False, output_format=False):
        GFAParser.save_graph(graph=self, output_path=output_file, force_format=output_format, minimal_graph=minimal)

