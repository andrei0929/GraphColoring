import random
import statistics
import sys
import time


def _generate_parent(length, geneSet, get_fitness):
    """

    Generate an array of random genes from geneSet of set length

    :param length: length of randomly generated genes array
    :param geneSet: set of possible genes
    :param get_fitness: fitness function
    :return: random genes array of length length
    """

    genes = []
    while len(genes) < length:
        sampleSize = min(length - len(genes), len(geneSet))
        genes.extend(random.sample(geneSet, sampleSize))
    fitness = get_fitness(genes)
    return Chromosome(genes, fitness)


def _mutate(parent, geneSet, get_fitness):
    """

    Mutates a random gene from the parent's genes

    :param parent: parent who's genes will be mutated
    :param geneSet: set of all possible genes
    :param get_fitness: fitness function
    :return: child with one mutated gene from parent
    """
    childGenes = parent.Genes[:]
    index = random.randrange(0, len(parent.Genes))
    newGene, alternate = random.sample(geneSet, 2)
    childGenes[index] = alternate if newGene == childGenes[index] else newGene
    fitness = get_fitness(childGenes)
    return Chromosome(childGenes, fitness)


def _mutate_custom(parent, custom_mutate, get_fitness):
    """

    Perform a custom mutation on the parent's genes

    :param parent: parent who's genes will be mutated
    :param custom_mutate: function which computes how the custom mutation is performed
    :param get_fitness: fitness function
    :return: child with custom mutated genes from parent
    """
    childGenes = parent.Genes[:]
    custom_mutate(childGenes)
    fitness = get_fitness(childGenes)
    return Chromosome(childGenes, fitness)


def get_best(get_fitness, targetLen, optimalFitness, geneSet, display,
             custom_mutate=None):
    """

    Function which computes the best child which is of optimalFitness
    using repeated generation and mutation of an initial guess

    :param get_fitness: fitness function
    :param targetLen: number of genes to use when creating a new gene sequence
    :param optimalFitness: optimal fitness value
    :param geneSet: set of possible genes
    :param display: display function
    :param custom_mutate: custom mutation function
    :return: best child with fitness value equal to optimalFitness
    """
    if custom_mutate is None:
        def fnMutate(parent):
            return _mutate(parent, geneSet, get_fitness)
    else:
        def fnMutate(parent):
            return _mutate_custom(parent, custom_mutate, get_fitness)

    def fnGenerateParent():
        return _generate_parent(targetLen, geneSet, get_fitness)

    for improvement in _get_improvement(fnMutate, fnGenerateParent):
        display(improvement)
        if not optimalFitness > improvement.Fitness:
            return improvement


def _get_improvement(new_child, generate_parent):
    """

    Generator function which successively computes a mutation of the best previous guess
    and retains the one with better fitness

    :param new_child: function which states how should a new child be computed
    :param generate_parent: function which generated the first initial guess
    :return: the guess with the larger fitness value
    """
    bestParent = generate_parent()
    yield bestParent
    while True:
        child = new_child(bestParent)
        if bestParent.Fitness > child.Fitness:
            continue
        if not child.Fitness > bestParent.Fitness:
            bestParent = child
            continue
        yield child
        bestParent = child


class Chromosome:
    def __init__(self, genes, fitness):
        self.Genes = genes
        self.Fitness = fitness
