class UnionFind:
    """Union-Find (disjoint set) data structure allowing to address by vertex name"""

    def __init__(self, vertex_names):
        #: Node name to id mapping
        self._name_to_id = {v: i for i, v in enumerate(vertex_names)}
        #: Pointer to the containing sets
        self._id = list(range(len(vertex_names)))
        #: Size of the set (_sz[_id[v]] is the size of the set that contains v)
        self._sz = [1] * len(vertex_names)
        #: Node id to name mapping
        self.id_to_name = {i: v for i, v in enumerate(vertex_names)}

    def name_of_id(self, i):
        return self.id_to_name[self.find(i)]

    def find(self, v):
        # assert type(v) is int
        j = v
        while j != self._id[j]:
            self._id[j] = self._id[self._id[j]]
            j = self._id[j]
        return j

    def find_by_name(self, v_name):
        return self.find(self._name_to_id[v_name])

    def union_by_name(self, v_name, w_name):
        self.union(self.find_by_name(v_name), self.find_by_name(w_name))

    def union(self, v, w):
        # assert type(v) is int
        # assert type(w) is int
        i = self.find(v)
        j = self.find(w)
        if i == j:
            return
        if self._sz[i] < self._sz[j]:
            self._id[i] = j
            self._sz[j] += self._sz[i]
        else:
            self._id[j] = i
        self._sz[i] += self._sz[j]

    def get_connected_components(
        self, IDs=None, allow_singletons: bool = False
    ) -> list[set[int]]:
        """computes lists of sets of connected components"""
        ccomponents = dict()
        for id in IDs if IDs else self._name_to_id.keys():
            union = self.name_of_id(self.find_by_name(id))
            if union not in ccomponents:
                ccomponents[union] = {id}
            else:
                ccomponents[union].add(id)
        return [l for l in list(ccomponents.values()) if len(l) > 1 or allow_singletons]
