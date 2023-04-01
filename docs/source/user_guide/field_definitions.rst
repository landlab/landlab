.. _standard_name_definitions:

======================
List of Landlab Fields
======================

The following tables list all of the field names used by *landlab* components. The *Used By* tab lists
component that use a given field (i.e. it is an input for the component), while the *Provided By* tab lists
components that provide the field as output.

.. jinja:: llcats

    .. tab:: Description

        .. list-table::
            :widths: 50 50
            :header-rows: 0

            {% for name, info in fields |dictsort %}
            * - *{{name}}*
              - {{info['desc']}}
            {% endfor %}

    .. tab:: Used By

        .. list-table::
            :widths: 50 50
            :header-rows: 0

            {% for name, info in fields |dictsort %}
            * - *{{name}}*
              - {% if info['used_by'] %}{% for cls in info['used_by'] %}:py:class:`~{{cls}}` {% endfor %}{% endif %}
            {% endfor %}

    .. tab:: Provided By

        .. list-table::
            :widths: 50 50
            :header-rows: 0

            {% for name, info in fields |dictsort %}
            * - *{{name}}*
              - {% if info['provided_by'] %}{% for cls in info['provided_by'] %}:py:class:`~{{cls}}` {% endfor %}{% endif %}
            {% endfor %}
