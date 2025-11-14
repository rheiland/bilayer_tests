# bilayer_tests

Create a sequence of simple bilayer models and test cell division and overall mechanics.

## Horizontal bilayer v0
* 3 cell types: basal, luminal, endpt
* only the luminal cells divide
* in custom_modules/custom.cpp, `custom_division_function` constrains where daughter cells are positioned
* the endpt cells exist to simplify the logic; they do not divide
* the "links" represent neighbors (not spring attachments). Rf View -> Plot options: Graph display: neighbors
<img src=.\images\horiz_bilayer_v0.gif width="50%">

