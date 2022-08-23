.. _standard_name_definitions:

Landlab Fields
==============

.. jinja:: llcats

    .. tab:: Description
    
        .. list-table::
            :widths: 50 50
            :header-rows: 1
           
            * - Field Name
              - Description
        
            {% for name, info in fields |dictsort %}
            * - {{name}}
              - {{info['desc']}}
            {% endfor %}
    
    .. tab:: Used By
    
        .. list-table::
            :widths: 50 50
            :header-rows: 1
            
            * - Field Name
              - Description
            
            {% for name, info in fields |dictsort %}
            * - {{name}}
              - {% if info['used_by'] %}{% for cls in info['used_by'] %}:py:class:`~{{cls}}` {% endfor %}{% endif %}
            {% endfor %}
            
    .. tab:: Provided By
    
        .. list-table::
            :widths: 50 50
            :header-rows: 1
            
            * - Field Name
              - Description
            
            {% for name, info in fields |dictsort %}
            * - {{name}}
              - {% if info['provided_by'] %}{% for cls in info['provided_by'] %}:py:class:`~{{cls}}` {% endfor %}{% endif %}
            {% endfor %}
    
