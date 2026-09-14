import json


class Structure:
    def __init__(self):
        self.nucleotides = {}
        self.tetrads = []
        self.tetrads_order = []
        self.single_tetrads = []
        self.tracts = []
        self.loops = []
        self.bulges = []
        self.strand_polarities = []
        self.paths = []
        self.basePairs = []

    def addNucleotide(self, name, data):
        self.nucleotides[name] = data

    def fromFile(self, path):
        with open(path) as file:
            data = json.load(file)
        return self.fromJsonDict(data)

    def fromString(self, json_string):
        data = json.loads(json_string)
        return self.fromJsonDict(data)

    def OrderTetrads(self, tetrad_pairs):
        tetrad_pairs = list(tetrad_pairs)
        tetrad_ordered = []

        if tetrad_pairs != None and len(tetrad_pairs) > 0:
            pair = tetrad_pairs.pop(0)
            tetrad_ordered.append(pair["tetrad1"])
            tetrad_ordered.append(pair["tetrad2"])
            to_remove = 0

            while len(tetrad_pairs) > 0:
                found = False
                for index_pair, pair in enumerate(tetrad_pairs):
                    for index, tetrad in enumerate(tetrad_ordered):
                        if tetrad == pair["tetrad1"]:
                            tetrad_ordered.insert(index + 1, pair["tetrad2"])
                            found = True
                            break
                        elif tetrad == pair["tetrad2"]:
                            tetrad_ordered.insert(index, pair["tetrad1"])
                            found = True
                            break
                    if found == True:
                        to_remove = index_pair
                        break
                tetrad_pairs.pop(to_remove)
                if len(tetrad_pairs) == 0:
                    break

        tetrad_ordered.reverse()
        return tetrad_ordered

    def fromJsonDict(self, json_dict):
        for data in json_dict["nucleotides"]:
            self.addNucleotide(data["fullName"], data)

        self.basePairs = json_dict.get("basePairs", [])

        for _, helice in enumerate(json_dict["helices"]):
            single_tetrads_local = []
            tetrad_unordered = {}
            tracts_all = []
            loops_all = []
            bulges_all = []
            strand_polarities_all = []
            paths_all = []

            for _, quadruplex in enumerate(helice["quadruplexes"]):
                tetrad_unordered_local = {}

                for data in quadruplex["tetrads"]:
                    tetrad_unordered[data["id"]] = data
                    tetrad_unordered_local[data["id"]] = data

                single_tetrads_local.append(tetrad_unordered_local)
                tracts_all.append(quadruplex.get("tracts", list()))
                loops_all.append(quadruplex.get("loops", list()))
                bulges_all.append(quadruplex.get("bulges", list()))
                strand_polarities_all.append(quadruplex.get("strandPolarities", list()))
                paths_all.append(quadruplex.get("path", list()))

            tetrad_ordered = self.OrderTetrads(helice["tetradPairs"])

            if len(tetrad_ordered) > 1:
                self.tetrads_order.append(tetrad_ordered)
                self.tetrads.append(tetrad_unordered)
            else:
                self.tetrads_order.append(list(tetrad_unordered.keys()))
                self.tetrads.append(tetrad_unordered)

            self.single_tetrads.append(single_tetrads_local)
            self.tracts.append(tracts_all)
            self.loops.append(loops_all)
            self.bulges.append(bulges_all)
            self.strand_polarities.append(strand_polarities_all)
            self.paths.append(paths_all)

        return self
