# bilayer_tests

Create a sequence of simple bilayer models and test cell division and overall mechanics.

## Horizontal bilayer v0
* 3 cell types: basal, luminal, endpt
* only the luminal cells divide
* in custom_modules/custom.cpp, `custom_division_function` constrains where daughter cells are positioned
* the endpt cells exist to simplify the logic; they do not divide
* the "links" represent neighbors (not spring attachments). Rf View -> Plot options: Graph display: neighbors
<img src=.\images\horiz_bilayer_v0.gif width="50%">

<img src=.\images\horiz_bilayer_v0_9days.png width="30%">
Letting it run longer shows less desirable results. Why? Probably a combination of not using spring
mechanics and allowing uninhibited, stochastic cell division. So we incrementally address each of
these in subsequent versions of the model.

