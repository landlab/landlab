# List of Components

````{jinja} llcats
```{{"{"}}list-table{{"}"}}
:widths: 50 50
:header-rows: 0

{% for name, component in components |dictsort %}
{% set parts = component['name'].split(".") %}
* - {{"{"}}class{{"}"}}`{{name}} <{{ component['name'] }}>`
  - {{ component['summary'] }}
{% endfor %}
```
````
