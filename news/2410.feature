Refactored `DataRecord.add_record` and `DataRecord.add_item` to avoid
`xarray.merge` when adding records or items. Missing item locations now use
different fill values for `element_id` and `grid_element`, preserving data
types.
