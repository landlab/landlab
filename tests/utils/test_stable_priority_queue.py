#! /usr/bin/env python

import numpy as np
import pytest

from landlab.utils import StablePriorityQueue


def test_add_subtract_examine():
    q = StablePriorityQueue()
    q.add_task("b", priority=2)
    q.add_task("a", priority=1)
    q.add_task(0, priority=0)
    q.add_task("c", priority=2)
    q.remove_task(0)
    assert q.pop_task() == "a"

    assert q.peek_at_task() == "b"

    assert np.all(q.tasks_currently_in_queue() == np.array(["b", "c"]))

    assert q.pop_task() == "b"

    assert np.all(q.tasks_ever_in_queue() == np.array(["b", "a", "0", "c"]))


def test_type_return_1():
    q = StablePriorityQueue()
    q.add_task(2, priority=2)
    q.add_task(1, priority=1)
    assert np.issubdtype(q.tasks_currently_in_queue().dtype, np.integer)


def test_type_return_2():
    q = StablePriorityQueue()
    q.add_task(np.pi)
    assert np.issubdtype(q.tasks_currently_in_queue().dtype, np.floating)


def test_empty_pop():
    q = StablePriorityQueue()
    with pytest.raises(KeyError):
        q.pop_task()


def test_empty_peek():
    q = StablePriorityQueue()
    with pytest.raises(KeyError):
        q.peek_at_task()


def test_inf_load_error():
    q = StablePriorityQueue()
    with pytest.raises(ValueError):
        q.add_task(float("inf"))


def test_overwrite():
    q = StablePriorityQueue()
    q.add_task(0, priority=5)
    q.add_task(1, priority=1)
    q.add_task(0, priority=0)
    assert q.pop_task() == 0
    assert len(q.tasks_currently_in_queue()) == 1
