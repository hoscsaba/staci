#!/usr/bin/env python

# Finding connected components in a bidirectional graph.
# By Mario Vilas (mvilas at gmail dot com)
# http://breakingcode.wordpress.com/2013/04/08/finding-connected-components-in-a-graph

# The graph nodes.
class Data(object):
    def __init__(self, name):
        self.__name  = name
        self.__links = set()

    @property
    def name(self):
        return self.__name

    @property
    def links(self):
        return set(self.__links)

    def add_link(self, other):
        self.__links.add(other)
        other.__links.add(self)

# The function to look for connected components.
def connected_components(nodes):

    # List of connected components found. The order is random.
    result = []

    # Make a copy of the set, so we can modify it.
    nodes = set(nodes)

    # Iterate while we still have nodes to process.
    while nodes:

        # Get a random node and remove it from the global set.
        n = nodes.pop()

        # This set will contain the next group of nodes connected to each other.
        group = {n}

        # Build a queue with this node in it.
        queue = [n]

        # Iterate the queue.
        # When it's empty, we finished visiting a group of connected nodes.
        while queue:

            # Consume the next item from the queue.
            n = queue.pop(0)

            # Fetch the neighbors.
            neighbors = n.links

            # Remove the neighbors we already visited.
            neighbors.difference_update(group)

            # Remove the remaining nodes from the global set.
            nodes.difference_update(neighbors)

            # Add them to the group of connected nodes.
            group.update(neighbors)

            # Add them to the queue, so we visit them in the next iterations.
            queue.extend(neighbors)

        # Add the group to the list of groups.
        result.append(group)

    # Return the list of groups.
    return result

# Load nodefile written by staci
def load_nodefile():

    a = Data("a")
    b = Data("b")
    a.add_link(b)   

    # Load node names
    fname = "nodelist.txt"
    with open(fname) as f:
        content = f.readlines()

    nodes=[]
    nodenames=[]
    i=1
    for elem in content:
        nodenames.append(elem.strip())
        nodes.append(Data(elem.strip()))
        i=i+1

    fname = "connected_nodes.txt"
    with open(fname) as f:
        content = f.readlines()

    for elem in content:
        node_pair = elem.split(",")
        if len(node_pair)>1:
                n1=node_pair[0].strip()
                n2=node_pair[1].strip()
                idx1=nodenames.index(n1)
                idx2=nodenames.index(n2)
                print "Node "+n1+" has index "+repr(idx1)
                print "Node "+n2+" has index "+repr(idx2)
                a=Data("a")
                b=Data("b")
                a.add_link(b)
                # nodes[idx1].add_link(nodes[idx2])
                
    return nodes

def connected_nodes():
    # # The first group, let's make a tree.
    a = Data("a")
    b = Data("b")
    # c = Data("c")
    # d = Data("d")
    # e = Data("e")
    # f = Data("f")
    a.add_link(b)    #      a
    # a.add_link(c)    #     / \
    # b.add_link(d)    #    b   c
    # c.add_link(e)    #   /   / \
    # c.add_link(f)    #  d   e   f
    # f.add_link(g)

    # # The second group, let's leave a single, isolated node.
    # g = Data("g")

    # # The third group, let's make a cycle.
    # h = Data("h")
    # i = Data("i")
    # j = Data("j")
    # k = Data("k")
    # h.add_link(i)    #    h----i
    # i.add_link(j)    #    |    |
    # j.add_link(k)    #    |    |
    # k.add_link(h)    #    k----j

    # Put all the nodes together in one big set.
    # nodes = {a, b, c, d, e, f, g, h, i, j, k}

    nodes = load_nodefile()

    # Find all the connected components.
    number = 1
    for components in connected_components(nodes):
        names = sorted(node.name for node in components)
        names = ", ".join(names)
        print "Group #%i: %s" % (number, names)
        number += 1

    # You should now see the following output:
    # Group #1: a, b, c, d, e, f
    # Group #2: g
    # Group #3: h, i, j, k

    # The test code...
if __name__ == "__main__":
    connected_nodes()