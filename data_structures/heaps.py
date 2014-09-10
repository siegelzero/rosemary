from collections import deque


################################################################################
# Leftist heaps
################################################################################


class _LeftistHeapNode(object):
    """
    Leftist heap node object.

    Each node has the following attributes:

        * value: The value of the node used for ordering.

        * rank: The distance from the node to the nearest leaf in the tree.
        Leftist trees are arranged so that the subtree with the shortest path to
        a leaf is on the right descendant.

        * left: Pointer to the left child.

        * right: Pointer to the right child.
    """
    __slots__ = ('value', 'rank', 'left', 'right')

    def __init__(self, value):
        self.value = value
        self.rank = 0
        self.left = None
        self.right = None


def _leftist_heap_merge(nodeA, nodeB):
    """
    Merges the two leftist heaps with given root nodes.

    Our implementation follows that given in Section 4.4.4 of "Graphs,
    Algorithms, and Optimization" by Kocay and Kreher. See also "Data Structures
    and Network Algorithms" by Tarjan for more detailed implementation notes.
    """
    if nodeA is None:
        return nodeB
    if nodeB is None:
        return nodeA

    # Make it so that nodeA has the smaller root.
    if nodeA.value > nodeB.value:
        nodeA, nodeB = nodeB, nodeA

    # If nodeA has no right child, we simply attach nodeB.
    if nodeA.right is None:
        nodeC = nodeB
    else:
        nodeC = _leftist_heap_merge(nodeA.right, nodeB)

    if nodeA.left is None:
        nodeA.right = None
        nodeA.left = nodeC
        nodeA.rank = 0
    else:
        # Make it so that the right child has the smaller rank.
        if nodeC.rank <= nodeA.left.rank:
            nodeA.right = nodeC
        else:
            nodeA.right = nodeA.left
            nodeA.left = nodeC

        nodeA.rank = nodeA.right.rank + 1

    return nodeA


class LeftistHeap(object):
    """
    Leftist heap object.

    A leftist heap is a binary tree with a priority queue ordering with a rank
    at every node. This rank gives the distance from the node to the nearest
    leaf. Leftist heaps are arranged so that the subtree with the shortest path
    to a leaf is on the right. The "leftist" property is that at each node, the
    depth of the right subtree never exceeds that of the left subtree.

    Leftist heaps have the same performance bounds as binary heaps for insert
    and delete_min, but allow merging to be performed in logarithmic time. For
    this reason, we typically only use leftist heaps when merges are needed.

    Our implementation follows that given in "Graphs, Algorithms, and
    Optimization" by Kocay and Kreher. For more details, see "Data Structures
    and Network Algorithms" by Tarjan and "Data Structures and Algorithm
    Analysis" by Weiss.
    """
    __slots__ = ('root', 'size')

    def __init__(self, values):
        """
        Creates a new leftist heap contining the given values.

        Input:
            * values: iterable (list, set, tuple)
                A list of values to insert into the heap.

        Details:
            Given a list of n values, we can create a leftist heap out of them
            in time O(n). We make each value into a one-item heap, and insert
            these into a queue. To create the heap, we remove the first two
            items from the queue, merge then, and add the resulting heap to the
            end of the queue. We repeat this until the queue contains only one
            heap, which we return.
        """
        num_values = len(values)
        size = 0
        nodes = deque([_LeftistHeapNode(val) for val in values])

        while num_values > 1:
            nodeA = nodes.popleft()
            nodeB = nodes.popleft()
            nodeA = _leftist_heap_merge(nodeA, nodeB)
            nodes.append(nodeA)
            num_values -= 1
            size += 1

        self.root = nodeA
        self.size = size + 1

    def __repr__(self):
        return "Leftist heap of size {}.".format(self.size)

    def __len__(self):
        return self.size

    def __iter__(self):
        size = self.size
        while size:
            value = self.delete_min()
            size -= 1
            yield value

    def insert(self, value):
        """
        Inserts the value into the heap, modifying self.

        Input:
            * value: value to insert

        Details:
            To insert an item into a leftist heap, we make the item a one-node
            heap and merge it with the existing heap. This requires O(log n)
            time.
        """
        nodeA = _LeftistHeapNode(value)
        self.size += 1

        if self.root is None:
            self.root = nodeA
        else:
            self.root = _leftist_heap_merge(self.root, nodeA)

    def find_min(self):
        """
        Returns the minimum value of the heap, without modifying self.

        Output:
            * value: minimum value in the heap.

        Details:
            The root of the leftist heap is always the smallest, so we can find
            this in constant time.
        """
        return self.root.value

    def delete_min(self):
        """
        Removes the minimum value from the heap and returns it, modifying self.

        Output:
            value: minimum value in the heap.

        Details:
            The minimum value of a leftist heap is the root of the tree. We
            remove this value and create the new heap by merging the left and
            right subtrees. This deletion takes O(log n) time, where n is the
            number of nodes in the two subtrees.
        """
        if self.root is None:
            value = None
        else:
            value = self.root.value
            self.root = _leftist_heap_merge(self.root.left, self.root.right)
            self.size -= 1
        return value

    def merge(self, other):
        """
        Merges the two leftist heaps, modifying self.

            * other: LeftistHeap object

        Details:
            Leftist heaps allow for efficient O(log n) merging of two trees.
            Suppose we want to merge two leftist trees A and B into a new tree
            T. We take the smaller of the roots of A and B to be the root of T.
            WLOG, assume that A is the smaller. Removing this element from A
            splits A into two subtrees, L and R. We need to make these three
            subtrees B, L, and R into two subtrees of T.  We do this by merging
            R and B into a new tree C, and then take C and L to be the left and
            right subtrees of T, placing the one with the smaller rank on the
            right.
        """
        self.root = _leftist_heap_merge(self.root, other.root)
        self.size += other.size


class _LazyLeftistHeapNode(object):
    __slots__ = ('value', 'rank', 'left', 'right', 'deleted')

    def __init__(self, value):
        self.value = value
        self.rank = 0
        self.left = None
        self.right = None
        self.deleted = False


def _purge(nodeA):
    if nodeA is None:
        return []

    if not nodeA.deleted:
        return [nodeA]
    else:
        return _purge(nodeA.left) + _purge(nodeA.right)


def _heapify(nodes):
    num_nodes = len(nodes)

    while num_nodes > 1:
        nodeA = nodes.popleft()
        nodeB = nodes.popleft()
        nodeA = _leftist_heap_merge(nodeA, nodeB)
        nodes.append(nodeA)
        num_nodes -= 1

    if nodes:
        return nodes.pop()
    else:
        return None


class _PairingHeapNode(object):
    def __init__(self, key, value=None):
        self.key = key
        self.value = value if value else key
        self.nextSibling = None
        self.leftChild = None
        self.prev = None


def _link_nodes(nodeA, nodeB):
    if nodeA is None:
        return nodeB
    if nodeB is None:
        return nodeA

    if nodeB.key < nodeA.key:
        nodeB.prev = nodeA.prev
        nodeA.prev = nodeB
        nodeA.nextSibling = nodeB.leftChild

        if nodeA.nextSibling is not None:
            nodeA.nextSibling.prev = nodeA

        nodeB.leftChild = nodeA

        return nodeB

    else:
        nodeB.prev = nodeA
        nodeA.nextSibling = nodeB.nextSibling

        if nodeA.nextSibling is not None:
            nodeA.nextSibling.prev = nodeA

        nodeB.nextSibling = nodeA.leftChild

        if nodeB.nextSibling is not None:
            nodeB.nextSibling.prev = nodeB

        nodeA.leftChild = nodeB

        return nodeA


def _combine_siblings(firstSibling):
    if firstSibling.nextSibling is None:
        return firstSibling

    node = firstSibling
    treeArray = []
    count = 0
    while node is not None:
        treeArray.append(node)
        count += 1
        node.prev.nextSibling = None
        node = node.nextSibling

    paired = []
    while len(treeArray) > 1:
        first = treeArray.pop()
        second = treeArray.pop()
        paired.append(_link_nodes(first, second))

    if treeArray:
        last = treeArray.pop()
    else:
        last = paired.pop()

    for tree in reversed(paired):
        last = _link_nodes(last, tree)

    return last


class PairingHeap(object):
    def __init__(self, node=None):
        self.root = node
        self.lookup = {}

    def insert(self, node):
        if self.root is None:
            self.root = node
        else:
            self.root = _link_nodes(self.root, node)

    def insert_values(self, values):
        for val in values:
            node = _PairingHeapNode(val)
            self.lookup[node.value] = node
            self.insert(node)

    def decrease_key(self, value, key):
        node = self.lookup[value]
        node.key = key

        if node != self.root:
            if node.nextSibling is not None:
                node.nextSibling.prev = node.prev

            if node.prev.leftChild == node:
                node.prev.leftChild = node.nextSibling
            else:
                node.prev.nextSibling = node.nextSibling

            node.nextSibling = None

            self.root = _link_nodes(self.root, node)

    def delete_min(self):
        root = self.root

        if root is None:
            return root

        if root.leftChild is None:
            self.root = None
        else:
            self.root = _combine_siblings(root.leftChild)

        del self.lookup[root.value]
        return root
