from collections import deque


class LBNode(object):
    __slots__ = ['value', 'rank', 'left', 'right']

    def __init__(self, value):
        self.value = value
        self.rank = 0
        self.left = None
        self.right = None


def _LBmerge(nodeA, nodeB):
    if nodeA is None:
        return nodeB
    if nodeB is None:
        return nodeA

    if nodeA.value > nodeB.value:
        nodeA, nodeB = nodeB, nodeA

    if nodeA.right is None:
        nodeC = nodeB
    else:
        nodeC = _LBmerge(nodeA.right, nodeB)

    if nodeA.left is None:
        nodeA.right = None
        nodeA.left = nodeC
        nodeA.rank = 0
    else:
        if nodeC.rank <= nodeA.left.rank:
            nodeA.right = nodeC
        else:
            nodeA.right = nodeA.left
            nodeA.left = nodeC

        nodeA.rank = nodeA.right.rank + 1

    return nodeA


class LBTree(object):
    __slots__ = 'root'

    def __init__(self):
        self.root = None

    def insert_values(self, values):
        nodes = deque([LBNode(val) for val in values])
        while len(nodes) > 1:
            nodeA = nodes.popleft()
            nodeB = nodes.popleft()
            nodeA = _LBmerge(nodeA, nodeB)
            nodes.append(nodeA)

        self.root = nodeA

    def insert(self, value):
        nodeA = LBNode(value)

        if self.root is None:
            self.root = nodeA
        else:
            self.root = _LBmerge(self.root, nodeA)

    def peek(self):
        return self.root.value

    def extract(self):
        if self.root is None:
            value = None
        else:
            value = self.root.value
            self.root = _LBmerge(self.root.left, self.root.right)
        return value

    def merge(self, other):
        self.root = _LBmerge(self.root, other.root)
