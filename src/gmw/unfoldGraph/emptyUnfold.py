from unfoldGraph.abstrctUnfold import AbstrctUnfolder
from mergeNodes.mergeBrother import BrotherMerger
from mergeNodes.splitParent import ParentSpliter
class EmptyUnfolder(AbstrctUnfolder):
    def unfold_graph(self):
        if self.merge_brother:
            merger = BrotherMerger(self.graph)
            merger.merge_brother()
            self._merge_nodes()
        if self.split_parent:
            spliter = ParentSpliter(self.graph)
            spliter.split_parent()
            self._merge_nodes()