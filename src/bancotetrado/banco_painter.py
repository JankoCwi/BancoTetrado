import math
import svgwrite


class Point:
    def __init__(self, x, y):

        self.x = x
        self.y = y

    def __iter__(self):

        yield self.x
        yield self.y

    def __str__(self):

        return "{0},{1}".format(self.x, self.y)

    def __repr__(self):

        return "Point({0}, {1})".format(self.x, self.y)


class Nucleotide:
    def __init__(self, name, data, row, column, polarity):

        self.name = name
        self.data = data
        self.row = row
        self.column = column
        self.polarity = polarity
        self.point = Point(0, 0)

        self.chain = data["chain"]
        self.index = data["index"]
        self.number = data["number"]
        self.short_name = data["shortName"]
        self.glycosidic_bond = data["glycosidicBond"]
        self.sugar_pucker = data.get("sugarPucker", "")



    def Letter(self):
        if self.glycosidic_bond == "syn":
            return self.short_name.lower()
        return self.short_name.upper()

    def FontWeight(self):
        if self.sugar_pucker == "South":
            return "bold"
        return "normal"

    def FontStyle(self):
        if self.sugar_pucker == "North":
            return "italic"
        return "normal"


class Diagram:
    def __init__(self, structure, helix_id, quadruplex_id):

        self.structure = structure
        self.helix_id = helix_id
        self.quadruplex_id = quadruplex_id
        self.tetrads = []
        self.tracts = structure.tracts[helix_id][quadruplex_id]
        self.loops = structure.loops[helix_id][quadruplex_id]
        self.bulges = structure.bulges[helix_id][quadruplex_id]
        self.strand_polarities = structure.strand_polarities[helix_id][quadruplex_id]
        self.path = structure.paths[helix_id][quadruplex_id]
        self.nucleotides = {}
        self.connections = []
        self.loop_labels = {}
        self.loop_inclusions = {}
        self.bulge_labels = {}
        self.starts = []
        self.ends = []

        self.PrepareTetrads()
        self.PrepareNucleotides()
        self.PrepareLoopLabels()
        self.PrepareBulgeLabels()
        self.PrepareConnections()
        self.PrepareEnds()

    def PrepareTetrads(self):

        tetrads = self.structure.single_tetrads[self.helix_id][self.quadruplex_id]
        order = self.structure.tetrads_order[self.helix_id]

        for name in order:
            if name in tetrads:
                self.tetrads.append(tetrads[name])



    def TractOrder(self):
        order = []

        for item in self.path:
            digits = ""
            for char in item:
                if char.isdigit():
                    digits += char
            if digits == "":
                continue

            tract_index = int(digits) - 1
            if tract_index >= 0 and tract_index < len(self.tracts) and tract_index not in order:
                order.append(tract_index)

        for index in range(len(self.tracts)):
            if index not in order:
                order.append(index)

        return order

    def PrepareNucleotides(self):
        column_by_name = {}
        polarity_by_name = {}

        for column, tract_index in enumerate(self.TractOrder()):
            tract = self.tracts[tract_index]
            for index, name in enumerate(tract):
                column_by_name[name] = column

                polarity = "plus"
                if len(self.strand_polarities) > tract_index and \
                   len(self.strand_polarities[tract_index]) > index:
                    polarity = self.strand_polarities[tract_index][index]
                polarity_by_name[name] = polarity

        for row, tetrad in enumerate(self.tetrads):

            for name in [tetrad["nt1"], tetrad["nt2"], tetrad["nt3"], tetrad["nt4"]]:
                
                data = self.structure.nucleotides[name]
                self.nucleotides[name] = Nucleotide(name, data, row, column_by_name[name], polarity_by_name.get(name, "plus"))

    def CoreNames(self):
        return set(self.nucleotides.keys())

    def OrderedCoreByChain(self):

        chains = {}

        for name, nucl in self.nucleotides.items():
            if nucl.chain not in chains:
                chains[nucl.chain] = []
            chains[nucl.chain].append(name)



        for chain in chains:
            chains[chain].sort(key = lambda name: self.nucleotides[name].index)

        return chains




    def FindCoreBeforeAfter(self, first_nucleotide, last_nucleotide):
       
        prev_name = ""
        next_name = ""
        prev_index = -1
        next_index = 999999999

        for name, nucl in self.nucleotides.items():
            if nucl.chain != first_nucleotide["chain"]:
                continue

            if nucl.index < first_nucleotide["index"] and nucl.index > prev_index:
                prev_index = nucl.index
                prev_name = name

            if nucl.index > last_nucleotide["index"] and nucl.index < next_index:
                next_index = nucl.index
                next_name = name

        return prev_name, next_name



    def LoopNucleotidePairsCore(self, name):
        core_names = self.CoreNames()

        for pair in self.structure.basePairs:
            if pair.get("inTetrad", False):
                continue

            if pair["nt1"] == name and pair["nt2"] in core_names:
                return True

            if pair["nt2"] == name and pair["nt1"] in core_names:
                return True

        return False


    def PrepareLoopLabels(self):

        for loop in self.loops:

            loop_nucleotides = loop["nucleotides"]
            first_loop = self.structure.nucleotides[loop_nucleotides[0]]
            last_loop = self.structure.nucleotides[loop_nucleotides[-1]]

            prev_name, next_name = self.FindCoreBeforeAfter(first_loop, last_loop)

            if prev_name != "" and next_name != "":

                if len(loop_nucleotides) == 1 and self.LoopNucleotidePairsCore(loop_nucleotides[0]):
                    self.loop_inclusions[(prev_name, next_name)] = loop_nucleotides
                else:
                    self.loop_labels[(prev_name, next_name)] = len(loop_nucleotides)

    def AddBulgeLabel(self, bulge_names):

        first_bulge = self.structure.nucleotides[bulge_names[0]]
        last_bulge = self.structure.nucleotides[bulge_names[-1]]

        prev_name, next_name = self.FindCoreBeforeAfter(first_bulge, last_bulge)

        if prev_name != "" and next_name != "":
            key = (prev_name, next_name)

            if key in self.loop_labels:
                return


            self.bulge_labels[(prev_name, next_name)] = bulge_names


    def PrepareBulgeLabels(self):
        core_names = self.CoreNames()


        for tract in self.tracts:
            for i in range(len(tract) - 1):

                name_a = tract[i]
                name_b = tract[i + 1]
                nucl_a = self.structure.nucleotides[name_a]
                nucl_b = self.structure.nucleotides[name_b]

                start_index = min(nucl_a["index"], nucl_b["index"])
                end_index = max(nucl_a["index"], nucl_b["index"])

                bulge_names = []

                for name, data in self.structure.nucleotides.items():
                    if data["chain"] == nucl_a["chain"] and \
                       data["index"] > start_index and \
                       data["index"] < end_index and \
                       name not in core_names:
                        bulge_names.append(name)

                if len(bulge_names) > 0:
                    bulge_names.sort(key = lambda name: self.structure.nucleotides[name]["index"])

                    if nucl_a["index"] < nucl_b["index"]:
                        key = (name_a, name_b)
                    else:
                        key = (name_b, name_a)


                    if len(bulge_names) == 1:
                        self.bulge_labels[key] = bulge_names
                    else:
                        self.loop_labels[key] = len(bulge_names)
                        self.bulge_labels[key] = bulge_names



        for bulge in self.bulges:
            if isinstance(bulge, str):
                self.AddBulgeLabel([bulge])

            elif "nucleotides" in bulge:
                self.AddBulgeLabel(bulge["nucleotides"])

    def PrepareConnections(self):
        chains = self.OrderedCoreByChain()

        for _, names in chains.items():
            for i in range(len(names) - 1):
                name_a = names[i]
                name_b = names[i + 1]
                label = 0
                label_type = ""
                bulge_names = []

                loop_inclusion_names = []

                if (name_a, name_b) in self.loop_labels:
                    label = self.loop_labels[(name_a, name_b)]
                    label_type = "loop"

                if (name_a, name_b) in self.loop_inclusions:
                    loop_inclusion_names = self.loop_inclusions[(name_a, name_b)]
                    if label_type == "":
                        label_type = "loop_inclusion"

                if (name_a, name_b) in self.bulge_labels:
                    label_type = "bulge"
                    bulge_names = self.bulge_labels[(name_a, name_b)]

                self.connections.append({
                    "from": name_a,
                    "to": name_b,
                    "label": label,
                    "label_type": label_type,
                    "bulge_names": bulge_names,
                    "loop_inclusion_names": loop_inclusion_names,
                })

    def PrepareEnds(self):
        chains = self.OrderedCoreByChain()

        for _, names in chains.items():
            if len(names) > 0:
                self.starts.append(names[0])
                self.ends.append(names[-1])


class BancoSvgMaker:
    def __init__(self, file_path, diagram):
        self.diagram = diagram
        self.file_path = file_path

        self.scale = 1.0
        self.font_family ="TimesNewRoman"

        self.margin = 72.0 * self.scale
        self.column_spacing = 82.0 * self.scale
        self.row_spacing = 72.0 * self.scale
        self.font_size = 34.0 * self.scale
        self.number_font_size = 9.0 * self.scale
        self.bulge_font_size = 22.0 * self.scale
        self.circle_font_size = 13.0 * self.scale
        self.circle_radius = 10.0 * self.scale
        self.stroke_width = 2.3 * self.scale
        self.nucleotide_radius = 16.0 * self.scale
        self.extra_nucleotide_radius = 9.0 * self.scale
        self.width = self.margin * 2.0 + self.column_spacing * 3.0
        self.height = self.margin * 2.0 + self.row_spacing * max(len(diagram.tetrads) - 1, 0)

        self.svg = svgwrite.Drawing(file_path, size = (self.width, self.height), profile = "full")

    def PositionNucleotides(self):
        for _, nucl in self.diagram.nucleotides.items():
            x = self.margin + nucl.column * self.column_spacing
            y = self.margin + nucl.row * self.row_spacing
            nucl.point = Point(x, y)

    def ConnectionPoint(self, point_a, point_b, t):
        return Point(point_a.x + (point_b.x - point_a.x) * t, \
                     point_a.y + (point_b.y - point_a.y) * t)

    def ShortenLine(self, point_a, point_b):
        return self.ShortenLineByRadii(point_a, point_b, self.nucleotide_radius, self.nucleotide_radius)

    def ShortenLineByRadii(self, point_a, point_b, radius_a, radius_b):
        dx = point_b.x - point_a.x
        dy = point_b.y - point_a.y
        length = math.sqrt(dx * dx + dy * dy)

        if length == 0:
            return point_a, point_b

        ux = dx / length
        uy = dy / length
        shift_a = min(radius_a, length * 0.30)
        shift_b = min(radius_b, length * 0.30)

        return Point(point_a.x + ux * shift_a, point_a.y + uy * shift_a), \
               Point(point_b.x - ux * shift_b, point_b.y - uy * shift_b)

    def DrawLine(self, point_a, point_b):
        line = self.svg.line(start = point_a, end = point_b, \
                stroke = "black", stroke_width = self.stroke_width)
        self.svg.add(line)

    def DrawCircleLabel(self, point, label):
        self.svg.add(self.svg.circle(point, r = self.circle_radius, \
                stroke = "black", stroke_width = self.stroke_width, \
                fill = "white"))

        text = self.svg.text(str(label), fill = "black", \
                insert = (point.x, point.y), \
                style = "text-anchor:middle;dominant-baseline:central;alignment-baseline:middle", \
                font_size = self.circle_font_size, \
                font_weight = "bold", \
                font_family = self.font_family)
        self.svg.add(text)

    def ExtraNucleotideLetter(self, data):
        if data["glycosidicBond"] == "syn":
            return data["shortName"].lower()
        return data["shortName"].upper()

    def ExtraNucleotideFontWeight(self, data):
        if data.get("sugarPucker", "") == "South":
            return "bold"
        return "normal"

    def ExtraNucleotideFontStyle(self, data):
        if data.get("sugarPucker", "") == "North":
            return "italic"
        return "normal"

    def DrawExtraNucleotide(self, name, point, font_size, upside_down = False):
        data = self.diagram.structure.nucleotides[name]
        transform = "translate({0}, {1})".format(point.x, point.y)
        if upside_down:
            transform += " rotate(180)"

        text = self.svg.text(self.ExtraNucleotideLetter(data), fill = "black", \
                transform = transform, \
                style = "text-anchor:middle;dominant-baseline:central;alignment-baseline:middle", \
                font_size = font_size, \
                font_weight = self.ExtraNucleotideFontWeight(data), \
                font_style = self.ExtraNucleotideFontStyle(data), \
                font_family = self.font_family)
        self.svg.add(text)

    def DrawBulge(self, connection, point_a, point_b):
        point = self.ConnectionPoint(point_a, point_b, 0.5)
        x_shift = 17.0 * self.scale
        y_shift = self.bulge_font_size * 0.5

        if connection["label"] > 0:
            symbol_x = point.x - self.circle_radius - x_shift * 0.4
        else:
            symbol_x = point.x - x_shift * 0.10

        symbol = self.svg.text("<", fill = "black", \
                insert = (symbol_x, point.y), \
                style = "text-anchor:middle;dominant-baseline:central;alignment-baseline:middle", \
                font_size = self.bulge_font_size, \
                font_weight = "bold", \
                font_family = self.font_family)
        self.svg.add(symbol)

        if connection["label"] > 0:
            return

        for index, name in enumerate(connection["bulge_names"]):
            y = point.y + (index - (len(connection["bulge_names"]) - 1) / 2.0) * self.bulge_font_size
            self.DrawExtraNucleotide(name, Point(point.x - x_shift, y + y_shift), self.bulge_font_size)


    def LoopInclusionPoints(self, connection, point_a, point_b):
        dx = point_b.x - point_a.x
        dy = point_b.y - point_a.y
        length = math.sqrt(dx * dx + dy * dy)
        if length == 0:
            return []

        nx = -dy / length
        ny = dx / length

        if ny < 0:
            nx = -nx
            ny = -ny

        offset = 18.0 * self.scale
        count = len(connection["loop_inclusion_names"])
        points = []

        for index, name in enumerate(connection["loop_inclusion_names"]):
            if count == 1:
                t = 0.35
            else:
                t = 0.30 + 0.18 * index


            base_point = self.ConnectionPoint(point_a, point_b, t)


            point = Point(base_point.x + nx * offset, base_point.y + ny * offset)


            for pair in self.diagram.structure.basePairs:

                if pair.get("inTetrad", False):
                    continue

                paired_name = ""

                if pair["nt1"] == name and pair["nt2"] in self.diagram.nucleotides:
                    paired_name = pair["nt2"]

                elif pair["nt2"] == name and pair["nt1"] in self.diagram.nucleotides:
                    paired_name = pair["nt1"]

                if paired_name != "":
                    point.y = self.diagram.nucleotides[paired_name].point.y

                    break




            points.append(point)

        return points

    def DrawLoopInclusions(self, connection, inclusion_points):
        for name, point in zip(connection["loop_inclusion_names"], inclusion_points):
            self.DrawExtraNucleotide(name, point, self.bulge_font_size)

    def DrawConnectionLabel(self, connection, point_a, point_b):

        if connection["label"] > 0:
            point = self.ConnectionPoint(point_a, point_b, 0.5)
            self.DrawCircleLabel(point, connection["label"])

        if connection["label_type"] == "bulge":
            self.DrawBulge(connection, point_a, point_b)


       

    def DrawBackground(self):
        self.svg.add(self.svg.rect(insert = (0, 0), size = (self.width, self.height), fill = "white"))

    def DrawConnections(self):
        for connection in self.diagram.connections:
            nucl_a = self.diagram.nucleotides[connection["from"]]
            nucl_b = self.diagram.nucleotides[connection["to"]]

            if len(connection["loop_inclusion_names"]) > 0:
                inclusion_points = self.LoopInclusionPoints(connection, nucl_a.point, nucl_b.point)
                route_points = [nucl_a.point] + inclusion_points + [nucl_b.point]

                for i in range(len(route_points) - 1):
                    radius_a = self.nucleotide_radius if i == 0 else self.extra_nucleotide_radius
                    radius_b = self.nucleotide_radius if i == len(route_points) - 2 else self.extra_nucleotide_radius
                    point_a, point_b = self.ShortenLineByRadii(route_points[i], route_points[i + 1], radius_a, radius_b)
                    self.DrawLine(point_a, point_b)

                self.DrawLoopInclusions(connection, inclusion_points)
                self.DrawConnectionLabel(connection, route_points[0], route_points[-1])
            else:
                point_a, point_b = self.ShortenLine(nucl_a.point, nucl_b.point)
                self.DrawLine(point_a, point_b)
                self.DrawConnectionLabel(connection, point_a, point_b)

    def DrawResidueNumber(self, nucl):
        pos = Point(nucl.point.x + self.font_size * 0.30, \
                    nucl.point.y + self.font_size * 0.35)
        text = self.svg.text(str(nucl.number), fill = "black", \
                insert = (pos.x, pos.y), \
                style = "text-anchor:start;dominant-baseline:central;alignment-baseline:middle", \
                font_size = self.number_font_size, \
                font_family = self.font_family)
        self.svg.add(text)

    def DrawNucleotide(self, nucl):
        transform = "translate({0}, {1})".format(nucl.point.x, nucl.point.y)
        if nucl.polarity == "minus":
            transform += " rotate(180)"

        text = self.svg.text(nucl.Letter(), fill = "black", \
                transform = transform, \
                style = "text-anchor:middle;dominant-baseline:middle", \
                font_size = self.font_size, \
                font_weight = nucl.FontWeight(), \
                font_style = nucl.FontStyle(), \
                font_family = self.font_family)
        self.svg.add(text)
        self.DrawResidueNumber(nucl)

    def EndDirection(self, name, label):
        for connection in self.diagram.connections:
            if label == "5'" and connection["from"] == name:
                point_a = self.diagram.nucleotides[connection["from"]].point
                point_b = self.diagram.nucleotides[connection["to"]].point
                return Point(point_a.x - point_b.x, point_a.y - point_b.y)

            if label == "3'" and connection["to"] == name:
                point_a = self.diagram.nucleotides[connection["from"]].point
                point_b = self.diagram.nucleotides[connection["to"]].point
                return Point(point_b.x - point_a.x, point_b.y - point_a.y)

        if label == "5'":
            return Point(-1.0, 0.0)
        return Point(1.0, 0.0)

    def DrawEndLabel(self, name, label):
        nucl = self.diagram.nucleotides[name]
        direction = self.EndDirection(name, label)
        length = math.sqrt(direction.x * direction.x + direction.y * direction.y)
        spacing = 32.0 * self.scale
        pos = Point(nucl.point.x + direction.x / length * spacing, \
                    nucl.point.y + direction.y / length * spacing)

        text = self.svg.text(label, fill = "black", \
                insert = (pos.x, pos.y), \
                style = "text-anchor:middle;dominant-baseline:central;alignment-baseline:middle", \
                font_size = self.font_size * 0.55, \
                font_weight = "bold", \
                font_family = self.font_family)
        self.svg.add(text)

    def DrawEnds(self):
        for name in self.diagram.starts:
            self.DrawEndLabel(name, "5'")
        for name in self.diagram.ends:
            self.DrawEndLabel(name, "3'")

    def DrawAll(self):
        self.PositionNucleotides()
        self.DrawBackground()
        self.DrawConnections()

        for _, nucl in self.diagram.nucleotides.items():
            self.DrawNucleotide(nucl)

        self.DrawEnds()


def Draw(structure, output_file):
    output_paths = []

    for helix_id in range(len(structure.single_tetrads)):
        quadruplex_count = len(structure.single_tetrads[helix_id])

        for quadruplex_id in range(quadruplex_count):
            diagram = Diagram(structure, helix_id, quadruplex_id)

            if quadruplex_count == 1:
                path = output_file + "_" + str(helix_id) + ".svg"
            else:
                path = output_file + "_" + str(helix_id) + "_" + str(quadruplex_id) + ".svg"

            svg_maker = BancoSvgMaker(path, diagram)
            svg_maker.DrawAll()
            svg_maker.svg.save(pretty = True)
            output_paths.append(path)
            print("Banco helix " + str(helix_id) + ", quadruplex " + str(quadruplex_id) + ": " + path)

    return output_paths
