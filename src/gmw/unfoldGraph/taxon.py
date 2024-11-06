class TaxonParser:
    def __init__(self, name_path, node_path):
        self.name_dict = {}
        self.node_dict = {}
        with open(name_path, "r") as file:
            for line in file:
                words = line.strip().split("\t|\t")
                name = words[1].lower().replace(" ", "_")
                if not name in self.name_dict:
                    self.name_dict[name] = words[0]
        with open(node_path, "r") as file:
            for line in file:
                words = line.strip().split("\t|\t")
                if not words[0] == words[1]:
                    self.node_dict[words[0]] = words[1]

    def _is_descendant(self, taxon, target):
        current = taxon
        while current in self.node_dict:
            if current == target:
                return True
            current = self.node_dict[current]
        return False

    def check_relationship(self, A, B):
        if self._is_descendant(A, B) or self._is_descendant(B, A):
            return True
        return False
    
    def name2taxid(self,name):
        lower_name = name.lower()
        if lower_name in self.name_dict:
            return self.name_dict[lower_name]
        raise ValueError(f"There is no scientific name \"{name}\"! Please correct the name or use tax_id instead.")