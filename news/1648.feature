
Added new, faster, aggregator functions for the `DataRecord`. These new functions,
`aggregate_items_at_link_count`, `aggregate_items_at_link_sum`, and
`aggregate_items_at_link_mean` are written in *cython* and are several orders
of magnitude faster than the previous method of using *merge*/*group_by*.
