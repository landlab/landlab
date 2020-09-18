#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 16:10:56 2020

@author: gtucker
"""

import numpy as np
from landlab.utils import StablePriorityQueueC


def test_stable_priority_queue_cython_implementation():
    q = StablePriorityQueueC()
    q.add_task('b', priority=2)
    q.add_task('a', priority=1)
    q.add_task(0, priority=0)
    q.add_task('c', priority=2)
    q.remove_task(0)
    task = q.peek_at_task()
    assert(task == 'b')
    task = q.pop_task()
    assert(task == 'a')
    task = q.peek_at_task()
    assert(task == 'b')
    assert(np.all(q.tasks_currently_in_queue() == np.array(['b', 'c'])))
    task = q.pop_task()
    assert(task == 'b')
