import random
import unittest
import datetime
import genetic


def load_data(localFileName):
    rules = set()
    nodes = set()
    with open(localFileName, mode='r') as infile:
        content = infile.read().splitlines()
    for row in content:
        if row[0] == 'e':  # e aa bb, aa and bb are node ids
            nodeIds = row.split(' ')[1:3]
            rules.add(Rule(nodeIds[0], nodeIds[1]))
            nodes.add(nodeIds[0])
            nodes.add(nodeIds[1])
            continue
        if row[0] == 'n':  # n aa ww, aa is a node id, ww is a weight
            nodeIds = row.split(' ')
            nodes.add(nodeIds[1])
    return rules, nodes


def get_fitness(genes, rules, stateIndexLookup):
    rulesThatPass = sum(1 for rule in rules
                        if rule.isValid(genes, stateIndexLookup))
    return rulesThatPass


def display(candidate, startTime):
    timeDiff = datetime.datetime.now() - startTime
    print("{0}\t{1}\t{2}".format(
        ''.join(map(str, candidate.Genes)),
        candidate.Fitness,
        str(timeDiff)))


def mutate(genes, geneset):
    # mutate between 1 and 4 genes
    count = random.randint(1, 4)
    while count > 0:
        indexA = random.randrange(0, len(genes))
        indexNewGeneA = random.randrange(0, len(geneset))
        genes[indexA] = geneset[indexNewGeneA]
        count -= 1


class GraphColoringTests(unittest.TestCase):

    def test_states(self):
        self.color("adjacent_states.col", ["Orange", "Yellow", "Green",
                                           "Blue"])

    def test_R100_1gb(self):
        self.color("R100_1gb.col",
                   ["Red", "Orange", "Yellow", "Green", "Blue", "Indigo"])

    def color(self, file, colors):
        rules, nodes = load_data(file)
        optimalValue = len(rules)
        colorLookup = {color[0]: color for color in colors}
        geneset = list(colorLookup.keys())
        startTime = datetime.datetime.now()
        nodeIndexLookup = {key: index
                           for index, key in enumerate(sorted(nodes))}

        def fnDisplay(candidate):
            display(candidate, startTime)

        def fnGetFitness(genes):
            return get_fitness(genes, rules, nodeIndexLookup)

        def fnMutate(genes):
            mutate(genes, geneset)

        best = genetic.get_best(fnGetFitness, len(nodes),
                                optimalValue, geneset, fnDisplay, custom_mutate=fnMutate)
        self.assertTrue(not optimalValue > best.Fitness)

        keys = sorted(nodes)
        for index in range(len(nodes)):
            print(keys[index] + " is " + colorLookup[best.Genes[index]])


class Rule:
    Node = None
    Adjacent = None

    def __init__(self, node, adjacent):
        if node < adjacent:
            node, adjacent = adjacent, node
        self.Node = node
        self.Adjacent = adjacent

    def __eq__(self, other):
        return self.Node == other.Node and \
               self.Adjacent == other.Adjacent

    def __hash__(self):
        return hash(self.Node) * 397 ^ hash(self.Adjacent)

    def __str__(self):
        return self.Node + " -> " + self.Adjacent

    # Verify if genes from 2 adjacent nodes are identical (colored with the same color)
    def isValid(self, genes, nodeIndexLookup):
        index = nodeIndexLookup[self.Node]
        adjacentStateIndex = nodeIndexLookup[self.Adjacent]

        return genes[index] != genes[adjacentStateIndex]


if __name__ == '__main__':
    unittest.main()
