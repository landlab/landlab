#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
cimport numpy as np
cimport cython

import heapq


cdef double LARGE_ELEV = 9999999999.0


cdef class StablePriorityQueueC:
    """Implements a stable priority queue, that tracks insertion order; i.e.,
    this is used to break ties.

    See https://docs.python.org/2/library/heapq.html#priority-queue-implementation-notes
    & https://www.sciencedirect.com/science/article/pii/S0098300413001337
    """
    cdef public list _pq
    cdef dict _entry_finder
    cdef float _REMOVED
    cdef int _counter

    def __init__(self):
        self._pq = []  # list of entries as a heap
        self._entry_finder = {}  # mapping of tasks to entries
        self._REMOVED = float("inf")  # placeholder for a removed task
        self._counter = 0  # unique sequence count

    def add_task(self, task, priority=0):
        """Add a new task or update the priority of an existing task."""
        if task == self._REMOVED:
            raise ValueError("StablePriorityQueue cannot accept tasks equal to INF!")
        if task in self._entry_finder:
            self.remove_task(task)
        #count = next(self._counter)
        self._counter += 1
        entry = [priority, self._counter, task]
        self._entry_finder[task] = entry
        heapq.heappush(self._pq, entry)

    def remove_task(self, task):
        """Mark an existing task as _REMOVED.

        Raise KeyError if not found.
        """
        entry = self._entry_finder.pop(task)
        entry[-1] = self._REMOVED

    def pop_task(self):
        """Remove and return the lowest priority task.

        Raise KeyError if empty.
        """
        while self._pq:
            priority, count, task = heapq.heappop(self._pq)
            #print("popped: ")
            #print(priority, count, task)
            if task != self._REMOVED:
                del self._entry_finder[task]
                return task
        raise KeyError("pop from an empty priority queue")

    def peek_at_task(self):
        """Return the lowest priority task without removal.

        Raise KeyError if empty.
        """
        counter = 0
        while self._pq:
            priority, count, task = self._pq[counter]
            #print('trying to peek at:')
            #print(priority, count, task)
            if task != self._REMOVED:
                return task
            counter += 1
        raise KeyError("peeked at an empty priority queue")

    def tasks_currently_in_queue(self):
        """Return array of nodes currently in the queue."""
        mynodes = [
            task for (priority, count, task) in self._pq if task != self._REMOVED
        ]
        return np.array(mynodes)
