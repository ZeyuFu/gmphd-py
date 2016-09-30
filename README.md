# GMPHD-py

This repository contains a Python implementation of the Gaussian Mixture - 
Probability Hypothesis Density (GM-PHD) filter described in [1] (forked from 
the [Dan Stowell repository](https://github.com/danstowell/gmphd)) and its 
application to underwater robotic mapping [2]. 

#### References

[1] B. N. Vo and W. K. Ma. *The Gaussian Mixture Probability Hypothesis Density Filter*
   IEEE Transactions on Signal Processing, 2006.

[2] T. Fabbri, F. Di Corato, D. Fenucci, D. Meucci and A. Caiti, 
*Multiple target tracking in seabed surveys using the GM-PHD filter* 
OCEANS 2015 - MTS/IEEE Washington, Washington, DC, 2015

[3] Stowell and M. D. Plumbley, *Multi-target pitch tracking of vibrato sources in
   noise using the GM-PHD filter*  In Proceedings of Proceedings of the 5th
   International Workshop on Machine Learning and Music, July 2012.

### Dependencies

GM-PHD Filter dependencies [3]:

- Numpy
- Scipy

Dependencies for the application of underwater robotic mapping [2]:

- MOOS
- MOOS-IVP
- Python-moos

### Notes

There are some differences from the GM-PHD algorithm described in Vo & Ma's paper:

- I have not implemented "spawning" of new targets from old ones, since I don't 
  need it. It would be straightforward to add it - see the original paper.

- Weights are adjusted at the end of pruning, so that pruning doesn't affect
  the total weight allocation.

- I provide an alternative approach to state-extraction (an alternative to
  Table 3 in the original paper) which makes use of the integral to decide how
  many states to extract.

### MOOSApp


### License

(C) 2016  Tommaso Fabbri - University of Pisa - Automation and Robotics Laboratory

This code represents free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
