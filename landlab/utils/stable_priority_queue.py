#!/usr/env/python

import heapq
import itertools
from collections import deque

import numpy as np


class StablePriorityQueue:
    """
    Implements a stable priority queue, that tracks insertion order; i.e., this
    is used to break ties.

    See https://docs.python.org/2/library/heapq.html#priority-queue-implementation-notes
    & https://www.sciencedirect.com/science/article/pii/S0098300413001337

    Examples
    --------
    >>> q = StablePriorityQueue()
    >>> q.add_task('b', priority=2)
    >>> q.add_task('a', priority=1)
    >>> q.add_task(0, priority=0)
    >>> q.add_task('c', priority=2)
    >>> q.remove_task(0)
    >>> q.pop_task()
    'a'
    >>> q.peek_at_task()
    'b'
    >>> np.all(q.tasks_currently_in_queue() == np.array(['b', 'c']))
    True
    >>> q.pop_task()
    'b'
    >>> np.all(q.tasks_ever_in_queue() == np.array(['b', 'a', '0', 'c']))
    True

    If only ints or floats are loaded into the array, the _in_queue methods
    return arrays with the corresponding data types:

    >>> q = StablePriorityQueue()
    >>> q.add_task(2, priority=2)
    >>> q.add_task(1, priority=1)
    >>> np.issubdtype(q.tasks_currently_in_queue().dtype, np.integer)
    True
    >>> q = StablePriorityQueue()
    >>> q.add_task(np.pi)
    >>> np.issubdtype(q.tasks_currently_in_queue().dtype, np.floating)
    True

    Popping from (or peeking at) an empty queue will throw a KeyError:

    >>> q = StablePriorityQueue()
    >>> try:
    ...     q.pop_task()
    ... except KeyError:
    ...     print('No tasks left')
    No tasks left
    """

    def __init__(self):
        self._pq = []  # list of entries as a heap
        self._entry_finder = {}  # mapping of tasks to entries
        self._REMOVED = float("inf")  # placeholder for a removed task
        self._counter = itertools.count()  # unique sequence count
        self._tasks_ever_in_queue = deque([])
        # last one tracks all nodes that have ever been added

    def add_task(self, task, priority=0):
        "Add a new task or update the priority of an existing task."
        if task == self._REMOVED:
            raise ValueError("StablePriorityQueue cannot accept tasks equal to INF!")
        if task in self._entry_finder:
            self.remove_task(task)
        count = next(self._counter)
        entry = [priority, count, task]
        self._entry_finder[task] = entry
        heapq.heappush(self._pq, entry)
        self._tasks_ever_in_queue.append(task)

    def remove_task(self, task):
        "Mark an existing task as _REMOVED.  Raise KeyError if not found."
        entry = self._entry_finder.pop(task)
        entry[-1] = self._REMOVED

    def pop_task(self):
        "Remove and return the lowest priority task. Raise KeyError if empty."
        while self._pq:
            priority, count, task = heapq.heappop(self._pq)
            if task is not self._REMOVED:
                del self._entry_finder[task]
                return task
        raise KeyError("pop from an empty priority queue")

    def peek_at_task(self):
        """
        Return the lowest priority task without removal. Raise KeyError if
        empty.
        """
        while self._pq:
            priority, count, task = self._pq[0]
            if task is not self._REMOVED:
                return task
        raise KeyError("peeked at an empty priority queue")

    def tasks_currently_in_queue(self):
        "Return array of nodes currently in the queue."
        mynodes = [
            task for (priority, count, task) in self._pq if task is not self._REMOVED
        ]
        return np.array(mynodes)

    def tasks_ever_in_queue(self):
        """
        Return array of all nodes ever added to this queue object. Repeats
        are permitted.
        """
        return np.array(self._tasks_ever_in_queue)
