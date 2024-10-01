Added new, faster, aggregator functions for the {class}`~.DataRecord`. These new functions,
{func}`~.aggregate_items_at_link_count`, {func}`~.aggregate_items_at_link_sum`, and
{func}`~.aggregate_items_at_link_mean` are written in *cython* and are several orders
of magnitude faster than the previous method of using *merge*/*group_by*.
