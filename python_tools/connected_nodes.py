class Data(object):
    def __init__(self, name):
        self.__name = name
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


def load_nodefile():

    a = Data("a")
    b = Data("b")
    a.add_link(b)

    # Load node names
    fname = "nodelist.txt"
    with open(fname) as f:
        content = f.readlines()

    nodes = []
    nodenames = []
    i = 1
    for elem in content:
        nodenames.append(elem.strip())
        nodes.append(Data(elem.strip()))
        i = i+1

    fname = "connected_nodes.txt"
    with open(fname) as f:
        content = f.readlines()

    for elem in content:
        node_pair = elem.split(",")
        if len(node_pair) > 1:
                n1 = node_pair[0].strip()
                n2 = node_pair[1].strip()
                idx1 = nodenames.index(n1)
                idx2 = nodenames.index(n2)
                nodes[idx1].add_link(nodes[idx2])
    return nodes


def connected_nodes():
    nodes = load_nodefile()

    # Find all the connected components.
    number = 1
    for components in connected_components(nodes):
        names = sorted(node.name for node in components)
        names = ", ".join(names)
        print("\nGroup #"+str(number)+":\n========= \n"+names)
        number += 1
    number -= 1

    if number > 1:
        print("\n\nAltogether "+str(number)+" groups.\n\n")
    else:
        print("\n\nAltogether "+str(number)+" group.\n\n") 

    # The test code...
if __name__ == "__main__":
    connected_nodes()
