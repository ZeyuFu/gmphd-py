# GMPHD-py

This is a Python implementation of the Gaussian Mixture - Probability Hypothesis Density (GM-PHD)
filter described in:

B. N. Vo and W. K. Ma. The Gaussian Mixture Probability Hypothesis Density Filter.
   IEEE Transactions on Signal Processing, 54(11):4091--4104, 2006.
   DOI: 10.1109/TSP.2006.881190

## Dependencies

- Numpy
- Scipy

The test file *gmphd-moos* requires python-moos
Tested with Python 2.7.


## Notes

There are some differences from the GM-PHD algorithm described in Vo & Ma's paper:

* I have not implemented "spawning" of new targets from old ones, since I don't 
  need it. It would be straightforward to add it - see the original paper.

* Weights are adjusted at the end of pruning, so that pruning doesn't affect
  the total weight allocation.

* I provide an alternative approach to state-extraction (an alternative to
  Table 3 in the original paper) which makes use of the integral to decide how
  many states to extract.


## License

gmphd-py is a python implementation of the Gaussian Mixture - Probability Hypothesis Density (GM-PHD) filter.
Copyright (C) 2016  Tommaso Fabbri - University of Pisa - Automation and Robotics Laboratory

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
